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

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.batch.item.ItemProcessor;

import uk.ac.ebi.eva.dbsnpimporter.exception.UndefinedHgvsAlleleException;
import uk.ac.ebi.eva.dbsnpimporter.models.SubSnpCoreFields;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Filters out all those variants with IUPAC ambiguity codes in the reference or alternate alleles.
 */
public class UnambiguousAllelesFilterProcessor implements ItemProcessor<SubSnpCoreFields, SubSnpCoreFields> {

    private static final Logger logger = LoggerFactory.getLogger(UnambiguousAllelesFilterProcessor.class);

    private Pattern pattern;

    public UnambiguousAllelesFilterProcessor() {
        pattern = Pattern.compile("[ACGT]*", Pattern.CASE_INSENSITIVE);
    }

    @Override
    public SubSnpCoreFields process(SubSnpCoreFields subSnpCoreFields) {
        String referenceInForwardStrand;
        try {
            referenceInForwardStrand = subSnpCoreFields.getReferenceInForwardStrand();
        } catch (UndefinedHgvsAlleleException hgvsReferenceUndefined) {
            logger.debug("Variant filtered out because reference allele is not defined: {} ({})", subSnpCoreFields,
                         hgvsReferenceUndefined);
            return null;
        }

        Matcher referenceMatcher = pattern.matcher(referenceInForwardStrand);
        if (!referenceMatcher.matches()) {
            logger.debug("Variant filtered out because reference allele is ambiguous: {}", subSnpCoreFields);
            return null;
        }

        String alternateInForwardStrand;
        try {
            alternateInForwardStrand = subSnpCoreFields.getAlternateInForwardStrand();
        } catch (UndefinedHgvsAlleleException hgvsAlternateUndefined) {
            logger.debug("Variant filtered out because alternate allele is not defined: {} ({})", subSnpCoreFields,
                         hgvsAlternateUndefined);
            return null;
        }

        Matcher alternateMatcher = pattern.matcher(alternateInForwardStrand);
        if (!alternateMatcher.matches()) {
            logger.debug("Variant filtered out because alternate allele is ambiguous: {}", subSnpCoreFields);
            return null;
        }

        return subSnpCoreFields;
    }

}
