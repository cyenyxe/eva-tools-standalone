#!/usr/bin/env python
import os
import re
import subprocess
import sys
from argparse import ArgumentParser, ArgumentError
from datetime import datetime

import yaml
from ebi_eva_common_pyutils import command_utils
from ebi_eva_common_pyutils.command_utils import run_command_with_output
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.config_utils import get_primary_mongo_creds_for_profile, get_accession_pg_creds_for_profile
from ebi_eva_common_pyutils.logger import logging_config, AppLogger
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query, execute_query

sys.path.append(os.path.dirname(__file__))

from remapping_config import load_config


class RemappingJob(AppLogger):

    @staticmethod
    def write_remapping_process_props_template(template_file_path):
        mongo_host, mongo_user, mongo_pass = get_primary_mongo_creds_for_profile(cfg['maven']['environment'],
                                                                                 cfg['maven']['settings_file'])
        pg_url, pg_user, pg_pass = get_accession_pg_creds_for_profile(cfg['maven']['environment'],
                                                                      cfg['maven']['settings_file'])
        with open(template_file_path, 'w') as open_file:
            open_file.write(f'''spring.datasource.driver-class-name=org.postgresql.Driver
spring.datasource.url={pg_url}
spring.datasource.username={pg_user}
spring.datasource.password={pg_pass}
spring.datasource.tomcat.max-active=3

spring.jpa.generate-ddl=true

spring.data.mongodb.host={mongo_host}
spring.data.mongodb.port=27017
spring.data.mongodb.database=eva_accession_sharded
spring.data.mongodb.username={mongo_user}
spring.data.mongodb.password={mongo_pass}

spring.data.mongodb.authentication-database=admin
mongodb.read-preference=secondaryPreferred
spring.main.web-environment=false
spring.main.allow-bean-definition-overriding=true
spring.jpa.properties.hibernate.jdbc.lob.non_contextual_creation=true
logging.level.uk.ac.ebi.eva.accession.remapping=INFO
parameters.chunkSize=1000
''')
        return template_file_path

    def get_job_information(self, assembly, taxid):
        query = (
            'SELECT source, scientific_name, assembly_accession, remapping_status, SUM(num_studies), '
            'SUM(num_ss_ids) '
            'FROM eva_progress_tracker.remapping_tracker '
            f"WHERE origin_assembly_accession='{assembly}' AND taxonomy='{taxid}' "
            'GROUP BY source, origin_assembly_accession, scientific_name, assembly_accession, remapping_status'
        )
        source_set = set()
        progress_set = set()
        scientific_name = None
        target_assembly = None
        n_study = 0
        n_variants = 0
        with get_metadata_connection_handle(cfg['maven']['environment'], cfg['maven']['settings_file']) as pg_conn:
            for source, scientific_name, target_assembly, progress_status, n_st, n_var in get_all_results_for_query(pg_conn, query):
                source_set.add(source)
                if progress_status:
                    progress_set.add(progress_status)
                n_study += n_st
                n_variants += n_var

        sources = ', '.join(source_set)
        if progress_set:
            progress_status = ', '.join(progress_set)
        else:
            progress_status = 'Pending'
        return sources, scientific_name, target_assembly, progress_status, n_study, n_variants

    def list_assemblies_to_process(self):
        query = 'SELECT DISTINCT origin_assembly_accession, taxonomy FROM eva_progress_tracker.remapping_tracker'
        with get_metadata_connection_handle(cfg['maven']['environment'], cfg['maven']['settings_file']) as pg_conn:
            for assembly, taxid in get_all_results_for_query(pg_conn, query):
                sources, scientific_name, target_assembly, progress_status, n_study, n_variants = \
                    self.get_job_information(assembly, taxid)
                print('\t'.join(str(e) for e in [sources, scientific_name, assembly, taxid, target_assembly, progress_status, n_study, n_variants]))

    def set_status_start(self, assembly, taxid):
        query = ('UPDATE eva_progress_tracker.remapping_tracker '
                 f"SET remapping_status='Started', remapping_start = '{datetime.now().isoformat()}' "
                 f"WHERE origin_assembly_accession='{assembly}' AND taxonomy='{taxid}'")
        with get_metadata_connection_handle(cfg['maven']['environment'], cfg['maven']['settings_file']) as pg_conn:
            execute_query(pg_conn, query)

    def set_status_end(self, assembly, taxid):
        query = ('UPDATE eva_progress_tracker.remapping_tracker '
                 f"SET remapping_status='Completed', remapping_end = '{datetime.now().isoformat()}' "
                 f"WHERE origin_assembly_accession='{assembly}' AND taxonomy='{taxid}'")
        with get_metadata_connection_handle(cfg['maven']['environment'], cfg['maven']['settings_file']) as pg_conn:
            execute_query(pg_conn, query)

    def set_status_failed(self, assembly, taxid):
        query = ('UPDATE eva_progress_tracker.remapping_tracker '
                 f"SET remapping_status = 'Failed' "
                 f"WHERE origin_assembly_accession='{assembly}' AND taxonomy='{taxid}'")
        with get_metadata_connection_handle(cfg['maven']['environment'], cfg['maven']['settings_file']) as pg_conn:
            execute_query(pg_conn, query)

    def set_counts(self, assembly, taxid, source, nb_variant_extracted=None, nb_variant_remapped=None,
                   nb_variant_ingested=None):
        set_statements = []
        query = (f"SELECT * FROM eva_progress_tracker.remapping_tracker "
                 f"WHERE origin_assembly_accession='{assembly}' AND taxonomy='{taxid}' AND source='{source}'")
        with get_metadata_connection_handle(cfg['maven']['environment'], cfg['maven']['settings_file']) as pg_conn:
            # Check that this row exists
            results = get_all_results_for_query(pg_conn, query)
        if results:
            if nb_variant_extracted is not None:
                set_statements.append(f"num_ss_extracted = {nb_variant_extracted}")
            if nb_variant_remapped is not None:
                set_statements.append(f"num_ss_remapped = {nb_variant_remapped}")
            if nb_variant_ingested is not None:
                set_statements.append(f"num_ss_ingested = {nb_variant_ingested}")

        if set_statements:
            query = ('UPDATE eva_progress_tracker.remapping_tracker '
                     'SET ' + ', '.join(set_statements) + ' '
                     f"WHERE origin_assembly_accession='{assembly}' AND taxonomy='{taxid}' AND source='{source}'")
            with get_metadata_connection_handle(cfg['maven']['environment'], cfg['maven']['settings_file']) as pg_conn:
                execute_query(pg_conn, query)

    def set_version(self, assembly, taxid, remapping_version=1):
        query = ('UPDATE eva_progress_tracker.remapping_tracker '
                 f"SET remapping_version='{remapping_version}' "
                 f"WHERE origin_assembly_accession='{assembly}' AND taxonomy='{taxid}'")

        with get_metadata_connection_handle(cfg['maven']['environment'], cfg['maven']['settings_file']) as pg_conn:
            execute_query(pg_conn, query)

    def check_processing_required(self, assembly, target_assembly, n_variants):
        if str(target_assembly) != 'None' and assembly != target_assembly and int(n_variants) > 0:
            return True
        return False

    def process_one_assembly(self, assembly, taxid, resume):
        self.set_status_start(assembly, taxid)
        base_directory = cfg['remapping']['base_directory']
        sources, scientific_name, target_assembly, progress_status, n_study, n_variants = self.get_job_information(assembly, taxid)
        if not self.check_processing_required(assembly, target_assembly, n_variants):
            self.info(f'Not Processing assembly {assembly} -> {target_assembly} for taxonomy {taxid}: '
                      f'{n_study} studies with {n_variants} '
                      f'found in {sources}')
            self.set_status_end(assembly, taxid)
            return

        self.info(f'Process assembly {assembly} for taxonomy {taxid}: {n_study} studies with {n_variants} '
                  f'found in {sources}')
        nextflow_remapping_process = os.path.join(os.path.dirname(__file__), 'remapping_process.nf')
        assembly_directory = os.path.join(base_directory, taxid, assembly)
        work_dir = os.path.join(assembly_directory, 'work')
        prop_template_file = os.path.join(assembly_directory, 'template.properties')
        os.makedirs(work_dir, exist_ok=True)
        remapping_log = os.path.join(assembly_directory, 'remapping_process.log')
        remapping_config_file = os.path.join(assembly_directory, 'remapping_process_config_file.yaml')
        remapping_config = {
            'taxonomy_id': taxid,
            'source_assembly_accession': assembly,
            'target_assembly_accession': target_assembly,
            'species_name': scientific_name,
            'output_dir': assembly_directory,
            'genome_assembly_dir': cfg['genome_downloader']['output_directory'],
            'template_properties': self.write_remapping_process_props_template(prop_template_file),
        }

        for part in ['executable', 'nextflow', 'jar']:
            remapping_config[part] = cfg[part]

        with open(remapping_config_file, 'w') as open_file:
            yaml.safe_dump(remapping_config, open_file)

        try:
            command = [
                cfg['executable']['nextflow'],
                '-log', remapping_log,
                'run', nextflow_remapping_process,
                '-params-file', remapping_config_file,
                '-work-dir', work_dir
            ]
            if resume:
                command.append('-resume')
            curr_working_dir = os.getcwd()
            os.chdir(assembly_directory)
            command_utils.run_command_with_output('Nextflow remapping process', ' '.join(command))
        except subprocess.CalledProcessError as e:
            self.error('Nextflow remapping pipeline failed')
            self.set_status_failed(assembly, taxid)
            raise e
        finally:
            os.chdir(curr_working_dir)
        self.set_status_end(assembly, taxid)
        self.count_variants_from_logs(assembly_directory, assembly, taxid)
        self.set_version(assembly, taxid)

    def count_variants_from_logs(self, assembly_directory, assembly, taxid):
        vcf_extractor_log = os.path.join(assembly_directory, 'logs', assembly + '_vcf_extractor.log')
        eva_remapping_count = os.path.join(assembly_directory, 'eva', assembly + '_eva_remapped_counts.yml')
        dbsnp_remapping_count = os.path.join(assembly_directory, 'dbsnp', assembly + '_dbsnp_remapped_counts.yml')
        eva_ingestion_log = os.path.join(assembly_directory, 'logs', assembly + '_eva_remapped.vcf_ingestion.log')
        dbsnp_ingestion_log = os.path.join(assembly_directory, 'logs', assembly + '_dbsnp_remapped.vcf_ingestion.log')

        eva_total, eva_written, dbsnp_total, dbsnp_written = count_variants_extracted(vcf_extractor_log)
        eva_candidate, eva_remapped, eva_unmapped = count_variants_remapped(eva_remapping_count)
        dbsnp_candidate, dbsnp_remapped, dbsnp_unmapped = count_variants_remapped(dbsnp_remapping_count)
        eva_ingested = count_variants_ingested(eva_ingestion_log)
        dbsnp_ingested = count_variants_ingested(dbsnp_ingestion_log)

        self.set_counts(
            assembly, taxid, 'EVA',
            nb_variant_extracted=eva_written,
            nb_variant_remapped=eva_remapped,
            nb_variant_ingested=eva_ingested
        )
        self.set_counts(
            assembly, taxid, 'DBSNP',
            nb_variant_extracted=dbsnp_written,
            nb_variant_remapped=dbsnp_remapped,
            nb_variant_ingested=dbsnp_ingested
        )

        self.info(f'For Taxonomy: {taxid} and Assembly: {assembly} Source: EVA ')
        self.info(f'Number of variant read:{eva_total}, written:{eva_written}, attempt remapping: {eva_candidate}, '
                  f'remapped: {eva_remapped}, failed remapped {eva_unmapped}')
        self.info(f'For Taxonomy: {taxid} and Assembly: {assembly} Source: DBSNP ')
        self.info(f'Number of variant read:{dbsnp_total}, written:{dbsnp_written}, attempt remapping: {dbsnp_candidate}, '
                  f'remapped: {dbsnp_remapped}, failed remapped {dbsnp_unmapped}')


def count_variants_remapped(count_yml_file):
    with open(count_yml_file) as open_file:
        data = yaml.safe_load(open_file)
    candidate_variants = data.get('all')
    remapped_variants = data.get('Flank_50', {}).get('Remapped', 0) + \
                        data.get('Flank_2000', {}).get('Remapped', 0) + \
                        data.get('Flank_50000', {}).get('Remapped', 0)
    unmapped_variants = data.get('filtered') or 0 + \
                        (data.get('Flank_50000', {}).get('total', 0) - data.get('Flank_50000', {}).get('Remapped', 0))
    return candidate_variants, remapped_variants, unmapped_variants


def parse_log_line(line):
    total = None
    written = None
    match = re.search(r'Items read = (\d+)', line)
    if match:
        total = int(match.group(1))
    match = re.search(r'items written = (\d+)', line)
    if match:
        written = int(match.group(1))
    return total, written


def count_variants_extracted(extraction_log):
    command = f'grep "EXPORT_EVA_SUBMITTED_VARIANTS_STEP" {extraction_log} | tail -1'
    log_line = run_command_with_output('Get total number of eva variants written', command, return_process_output=True)
    eva_total, eva_written = parse_log_line(log_line)
    command = f'grep "EXPORT_DBSNP_SUBMITTED_VARIANTS_STEP" {extraction_log} | tail -1'
    log_line = run_command_with_output('Get total number of dbsnp variants written', command, return_process_output=True)
    dbsnp_total, dbnp_written = parse_log_line(log_line)
    return eva_total, eva_written, dbsnp_total, dbnp_written


def count_variants_ingested(ingestion_log):
    command = f'grep "INGEST_REMAPPED_VARIANTS_FROM_VCF_STEP" {ingestion_log} | tail -1'
    log_line = run_command_with_output('Get total number of variants written', command, return_process_output=True)
    _, written = parse_log_line(log_line)
    return written


def main():
    argparse = ArgumentParser(description='Run entire variant remapping pipeline for a given assembly and taxonomy.')
    argparse.add_argument('--assembly', help='Assembly to be process')
    argparse.add_argument('--taxonomy_id', help='Taxonomy id to be process')
    argparse.add_argument('--list_jobs', help='Display the list of jobs to be run.', action='store_true', default=False)
    argparse.add_argument('--resume', help='If a process has been run already this will resume it.',
                          action='store_true', default=False)

    args = argparse.parse_args()

    load_config()

    if args.list_jobs:
        RemappingJob().list_assemblies_to_process()
    elif args.assembly and args.taxonomy_id:
        logging_config.add_stdout_handler()
        RemappingJob().process_one_assembly(args.assembly, args.taxonomy_id, args.resume)
    else:
        raise ArgumentError('One of (--assembly and --taxonomy_id) or --list_jobs options is required')


if __name__ == "__main__":
    main()
