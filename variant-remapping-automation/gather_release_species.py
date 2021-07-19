#!/usr/bin/env python
import json
import operator
import os
import re
from argparse import ArgumentParser
from collections import defaultdict

import atexit

import pandas as pd
import psycopg2.extras
import requests
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_common_pyutils.pg_utils import execute_query, get_all_results_for_query

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()


eutils_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
esearch_url = eutils_url + 'esearch.fcgi'
esummary_url = eutils_url + 'esummary.fcgi'
efetch_url = eutils_url + 'efetch.fcgi'
ensembl_url = 'http://rest.ensembl.org/info/assembly'


cache_file = 'cache.json'


def load_cache():
    global cache
    if os.path.exists(cache_file):
        with open(cache_file) as open_file:
            cache.update(json.load(open_file))


def save_cache():
    with open(cache_file, 'w') as open_file:
        json.dump(cache, open_file)


atexit.register(save_cache)


cache = defaultdict(dict)
load_cache()


def retrieve_assembly_summary_from_species_name(species):
    """Search for all ids of assemblies associated with a species by depaginating the results of the search query"""
    payload = {'db': 'Assembly', 'term': '"{}[ORGN]"'.format(species), 'retmode': 'JSON', 'retmax': 100}
    response = requests.get(esearch_url, params=payload)
    data = response.json()
    search_results = data.get('esearchresult', {})
    id_list = search_results.get('idlist', [])
    while int(search_results.get('retstart')) + int(search_results.get('retmax')) < int(search_results.get('count')):
        payload['retstart'] = int(search_results.get('retstart')) + int(search_results.get('retmax'))
        response = requests.get(esearch_url, params=payload)
        data = response.json()
        search_results = data.get('esearchresult', {})
        id_list += search_results.get('idlist', [])
    response = requests.get(esummary_url, params={'db': 'Assembly', 'id': ','.join(id_list), 'retmode': 'JSON'})
    summary_list = response.json()
    if summary_list and 'result' in summary_list:
        return [summary_list.get('result').get(uid) for uid in summary_list.get('result').get('uids')]


def most_recent_assembly(assembly_list):
    """Based on assembly summaries find the one submitted the most recently"""
    if assembly_list:
        return sorted(assembly_list, key=operator.itemgetter('submissiondate'))[-1]


def best_assembly(assembly_list):
    """Based on assembly summaries find the one with the highest scaffold N50"""
    if assembly_list:
        return sorted(assembly_list, key=operator.itemgetter('scaffoldn50'))[-1]


def retrieve_species_names_from_tax_id(taxid):
    """Search for a species scientific name based on the taxonomy id"""
    if str(taxid) not in cache['taxid_to_name']:
        logger.info(f'Query NCBI for taxonomy {taxid}', )
        payload = {'db': 'Taxonomy', 'id': taxid}
        r = requests.get(efetch_url, params=payload)
        match = re.search('<Rank>(.+?)</Rank>', r.text, re.MULTILINE)
        rank = None
        if match:
            rank = match.group(1)
        if rank not in ['species', 'subspecies']:
            logger.warning('Taxonomy id %s does not point to a species', taxid)
        match = re.search('<ScientificName>(.+?)</ScientificName>', r.text, re.MULTILINE)
        if match:
            cache['taxid_to_name'][str(taxid)] = match.group(1)
        else:
            logger.warning('No species found for %s' % taxid)
            cache['taxid_to_name'][str(taxid)] = None
    return taxid, cache['taxid_to_name'].get(str(taxid))


def retrieve_species_name_from_assembly_accession(assembly_accession):
    """Search for a species scientific name based on an assembly accession"""
    if assembly_accession not in cache['assembly_to_species']:
        logger.info(f'Query NCBI for assembly {assembly_accession}', )
        payload = {'db': 'Assembly', 'term': '"{}"'.format(assembly_accession), 'retmode': 'JSON'}
        data = requests.get(esearch_url, params=payload).json()
        if data and 'esearchresult' in data:
            assembly_id_list = data.get('esearchresult', {}).get('idlist', [])
            payload = {'db': 'Assembly', 'id': ','.join(assembly_id_list), 'retmode': 'JSON'}
            summary_list = requests.get(esummary_url, params=payload).json()
            all_species_names = set()
            for assembly_id in summary_list.get('result', {}).get('uids', []):
                assembly_info = summary_list.get('result').get(assembly_id)
                all_species_names.add((assembly_info.get('speciestaxid'), assembly_info.get('speciesname')))
            if len(all_species_names) == 1:
                cache['assembly_to_species'][assembly_accession] = all_species_names.pop()
            else:
                logger.warning('%s taxons found for assembly %s ' % (len(all_species_names), assembly_accession))
    return cache['assembly_to_species'].get(assembly_accession) or (None, None)


def retrieve_current_ensembl_assemblies(taxid_or_assembly):
    """
    Retrieve the assembly accession currently supported by ensembl for the provided taxid or assembly accession
    In both case it looks up the associated species name in NCBI and using the species name returns the currently
    supported assembly for this species.
    """
    logger.debug('Search for species name for %s', taxid_or_assembly)
    scientific_name = None
    if taxid_or_assembly and str(taxid_or_assembly).isdigit():
        # assume it is a taid
        taxid, scientific_name = retrieve_species_names_from_tax_id(taxid_or_assembly)
    elif taxid_or_assembly:
        # assume it is an assembly accession
        taxid, scientific_name = retrieve_species_name_from_assembly_accession(taxid_or_assembly)
    if scientific_name:
        logger.debug('Found %s', scientific_name)
        if scientific_name not in cache['scientific_name_to_ensembl']:
            logger.info(f'Query Ensembl for species {scientific_name}', )
            url = ensembl_url + '/' + scientific_name.lower().replace(' ', '_')
            response = requests.get(url, params={'content-type': 'application/json'})
            data = response.json()
            assembly_accession = str(data.get('assembly_accession'))
            cache['scientific_name_to_ensembl'][scientific_name] = assembly_accession
        return [str(taxid), str(scientific_name), cache['scientific_name_to_ensembl'].get(scientific_name)]

    return ['NA', 'NA', 'NA']


def find_all_eva_studies(accession_counts, private_config_xml_file):
    with get_metadata_connection_handle("development", private_config_xml_file) as pg_conn:
        query = (
            'SELECT DISTINCT a.vcf_reference_accession, pt.taxonomy_id, p.project_accession '
            'FROM project p '
            'LEFT OUTER JOIN project_analysis pa ON p.project_accession=pa.project_accession '
            'LEFT OUTER JOIN analysis a ON pa.analysis_accession=a.analysis_accession '
            'LEFT OUTER JOIN project_taxonomy pt ON p.project_accession=pt.project_accession '
            'WHERE p.ena_status=4 '   # Ensure that the project is public
            'ORDER BY pt.taxonomy_id, a.vcf_reference_accession'
        )
        data = []
        for assembly, tax_id, study in filter_studies(get_all_results_for_query(pg_conn, query)):
            taxid_from_ensembl, scientific_name, ensembl_assembly_from_taxid = retrieve_current_ensembl_assemblies(tax_id)
            _, _, ensembl_assembly_from_assembly = retrieve_current_ensembl_assemblies(assembly)

            count_ssid = 0
            if study in accession_counts:
                assembly_from_mongo, taxid_from_mongo, project_accession, count_ssid = accession_counts.pop(study)
                if assembly_from_mongo != assembly:
                    logger.error(
                        'For study %s, assembly from accessioning (%s) is different'
                        ' from assembly from metadata (%s) database.', study, assembly_from_mongo, assembly
                    )
                if taxid_from_mongo != tax_id:
                    logger.error(
                        'For study %s, taxonomy from accessioning (%s) is different'
                        ' from taxonomy from metadata (%s) database.', study, taxid_from_mongo, tax_id
                    )
            data.append({
                'Source': 'EVA',
                'Assembly': assembly,
                'Taxid': tax_id,
                'Scientific Name': scientific_name,
                'Study': study,
                'Number Of Variants (submitted variants)': count_ssid or 0,
                'Ensembl assembly from taxid': ensembl_assembly_from_taxid,
                'Ensembl assembly from assembly': ensembl_assembly_from_assembly,
                'Target Assembly': ensembl_assembly_from_taxid or ensembl_assembly_from_assembly
            })
    if len(accession_counts) > 0:
        logger.error('Accessioning database has studies (%s) absent from the metadata database', ', '.join(accession_counts))
    df = pd.DataFrame(data)
    df = df.groupby(
        ['Source', 'Assembly', 'Taxid', 'Scientific Name', 'Ensembl assembly from taxid',
         'Ensembl assembly from assembly', 'Target Assembly']
    ).agg(
        {'Study': 'count', 'Number Of Variants (submitted variants)': 'sum'}
    )
    df.rename(columns={'Study': 'number Of Studies'}, inplace=True)
    return df.reset_index()


def parse_accession_counts(accession_counts_file):
    accession_count = {}
    with open(accession_counts_file) as open_file:
        for line in open_file:
            sp_line = line.strip().split()
            accession_count[sp_line[0]] = int(sp_line[1])
    return accession_count


def get_accession_counts_per_study(private_config_xml_file, source):
    accession_count = {}
    with get_metadata_connection_handle("development", private_config_xml_file) as pg_conn:
        query = (
            'SELECT assembly_accession, taxid, project_accession, SUM(number_submitted_variants) '
            'FROM eva_stats.submitted_variants_load_counts '
            "WHERE source='%s'"
            'GROUP BY assembly_accession, taxid, project_accession ' % source
        )
        for assembly_accession, taxid, project_accession, count_ssid in get_all_results_for_query(pg_conn, query):
            accession_count[project_accession] = (assembly_accession, taxid, project_accession, count_ssid)
    return accession_count


def get_accession_counts_per_assembly(private_config_xml_file, source):
    accession_count = {}
    with get_metadata_connection_handle("development", private_config_xml_file) as pg_conn:
        query = (
            'SELECT assembly_accession, taxid, SUM(number_submitted_variants) '
            'FROM eva_stats.submitted_variants_load_counts '
            "WHERE source='%s'"
            'GROUP BY assembly_accession, taxid ' % source
        )
        for assembly_accession, taxid, count_ssid in get_all_results_for_query(pg_conn, query):
            accession_count[assembly_accession] = (assembly_accession, taxid, count_ssid)
    return accession_count


def filter_studies(query_results):
    """
    Remove studies from the EVA list that are either missing information (assembly or taxid)
    or that are human and therefore cannot be released
    """
    for assembly, tax_id, study in query_results:
        if not assembly or not tax_id:
            logger.error('Study %s is missing assembly (%s) or taxonomy id (%s)', study, assembly, tax_id)
        elif tax_id == 9606:
            logger.debug("Study %s is human and won't be released", study)
        else:
            yield assembly, tax_id, study


def parse_dbsnp_csv(input_file, accession_counts):
    """Parse the CSV file generated in the past year to get the DBSNP data"""
    df = pd.read_csv(input_file)
    df = df[df.Source != 'EVA']
    taxids = []
    scientific_names = []
    ensembl_assemblies_from_taxid = []
    ensembl_assemblies_from_assembly = []
    target_assemblies = []

    for index, record in df.iterrows():
        taxid, scientific_name, ensembl_assembly_from_taxid = retrieve_current_ensembl_assemblies(record['Taxid'])
        _, _, ensembl_assembly_from_assembly = retrieve_current_ensembl_assemblies(record['Assembly'])
        taxids.append(taxid)
        scientific_names.append(scientific_name)
        ensembl_assemblies_from_taxid.append(ensembl_assembly_from_taxid)
        ensembl_assemblies_from_assembly.append(ensembl_assembly_from_assembly)
        target_assemblies.append(ensembl_assembly_from_taxid or ensembl_assembly_from_assembly)
        if record['Assembly'] != 'Unmapped':
            _, _, count = accession_counts[record['Assembly']]
            if count != int(record['Number Of Variants (submitted variants)'].replace(',', '')):
                logger.error(
                    'Count in spreadsheet (%s) and in database (%s) are different for accession %s',
                    record['Number Of Variants (submitted variants)'], count, record['Assembly']
                )
    df['Scientific Name'] = scientific_names
    df['Ensembl assembly from taxid'] = ensembl_assemblies_from_taxid
    df['Ensembl assembly from assembly'] = ensembl_assemblies_from_assembly
    df['Target Assembly'] = target_assemblies
    df.replace(',', '', regex=True, inplace=True)
    return df


def create_table_for_progress(private_config_xml_file):
    with get_metadata_connection_handle("development", private_config_xml_file) as metadata_connection_handle:
        query_create_table = (
            'CREATE TABLE IF NOT EXISTS eva_progress_tracker.remapping_tracker '
            '(source TEXT, taxonomy INTEGER, scientific_name TEXT, origin_assembly_accession TEXT, num_studies INTEGER NOT NULL,'
            'num_ss_ids BIGINT NOT NULL, release_version INTEGER, assembly_accession TEXT, '
            'remapping_report_time TIMESTAMP DEFAULT NOW(), remapping_status TEXT, remapping_start TIMESTAMP, '
            'remapping_end TIMESTAMP, remapping_version TEXT, num_ss_extracted INTEGER, '
            'num_ss_remapped INTEGER, num_ss_ingested INTEGER, '
            'primary key(source, taxonomy, origin_assembly_accession, release_version))'
        )
    execute_query(metadata_connection_handle, query_create_table)


def insert_remapping_progress_to_db(private_config_xml_file, dataframe):
    list_to_remap = dataframe.values.tolist()
    if len(list_to_remap) > 0:
        with get_metadata_connection_handle("development", private_config_xml_file) as metadata_connection_handle:
            with metadata_connection_handle.cursor() as cursor:
                query_insert = (
                    'INSERT INTO eva_progress_tracker.remapping_tracker '
                    '(source, taxonomy, scientific_name, origin_assembly_accession, num_studies, '
                    'num_ss_ids, assembly_accession, release_version) '
                    'VALUES %s'
                )
                psycopg2.extras.execute_values(cursor, query_insert, list_to_remap)


def main():
    argparse = ArgumentParser(
        description='Gather the current set of studies from both EVA and dbSNP that can be remapped, clustered '
                    'and released. The source of EVA studies is the metadata database and the source of dbSNP studies '
                    "is last year's spreadsheet. The number of variants are populated from counts retrieved from "
                    ' eva_stats.')
    argparse.add_argument('--input', help='Path to the file containing the taxonomies and assemblies', required=True)
    argparse.add_argument('--output', help='Path to the file that will contain the input plus annotation',
                          required=True)
    argparse.add_argument('--private_config_xml_file', required=True,
                          help='Path to the file containing the username/passwords tp access '
                               'production and development databases')
    args = argparse.parse_args()
    output_header = ['Source', 'Taxid', 'Scientific Name', 'Assembly', 'number Of Studies',
                     'Number Of Variants (submitted variants)', 'Ensembl assembly from taxid',
                     'Ensembl assembly from assembly', 'Target Assembly']

    accession_counts_dbsnp = get_accession_counts_per_assembly(args.private_config_xml_file, 'dbSNP')
    df1 = parse_dbsnp_csv(args.input, accession_counts_dbsnp)
    accession_counts_eva = get_accession_counts_per_study(args.private_config_xml_file, 'EVA')
    df2 = find_all_eva_studies(accession_counts_eva, args.private_config_xml_file)
    df = pd.concat([df1, df2])
    df = df[output_header]
    df.to_csv(args.output, quoting=False, sep='\t', index=False)
    create_table_for_progress(args.private_config_xml_file)
    output_header = ['Source', 'Taxid', 'Scientific Name', 'Assembly', 'number Of Studies',
                     'Number Of Variants (submitted variants)', 'Target Assembly']
    df = df[output_header]
    df['Release'] = 3
    df = df[df['Source'] != 'DBSNP - filesystem']
    insert_remapping_progress_to_db(args.private_config_xml_file, df)


if __name__ == "__main__":
    # This script was copied from EVA2406: https://github.com/EBIvariation/eva-tasks/tree/master/tasks/eva_2406
    # It should be updated before being used beyond release 3
    main()
