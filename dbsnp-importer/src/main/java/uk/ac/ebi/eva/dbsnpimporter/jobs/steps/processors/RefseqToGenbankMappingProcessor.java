/*
 * Copyright 2017 EMBL - European Bioinformatics Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package uk.ac.ebi.eva.dbsnpimporter.jobs.steps.processors;

import org.springframework.batch.item.ItemProcessor;

import uk.ac.ebi.eva.commons.core.models.Region;
import uk.ac.ebi.eva.dbsnpimporter.contig.ContigMapping;
import uk.ac.ebi.eva.dbsnpimporter.models.SubSnpCoreFields;

public class RefseqToGenbankMappingProcessor implements ItemProcessor<SubSnpCoreFields, SubSnpCoreFields> {

    private ContigMapping contigMapping;

    public RefseqToGenbankMappingProcessor(ContigMapping contigMapping) {
        this.contigMapping = contigMapping;
    }

    @Override
    public SubSnpCoreFields process(SubSnpCoreFields subSnpCoreFields) throws Exception {
        if (!subSnpCoreFields.isValidRegion(subSnpCoreFields.getChromosomeRegion())
                && subSnpCoreFields.isValidRegion(subSnpCoreFields.getContigRegion())) {
            replaceContig(subSnpCoreFields);
        }
        return subSnpCoreFields;
    }

    private void replaceContig(SubSnpCoreFields subSnpCoreFields) {
        Region contigRegion = subSnpCoreFields.getContigRegion();
        String genbank = contigMapping.getGenbankOrDefault(contigRegion.getChromosome());
        contigRegion.setChromosome(genbank);
        subSnpCoreFields.setContigRegion(contigRegion);
    }
}
