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

import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import uk.ac.ebi.eva.commons.core.models.IVariantSource;
import uk.ac.ebi.eva.commons.core.models.pedigree.Sex;
import uk.ac.ebi.eva.dbsnpimporter.models.DbsnpBatch;
import uk.ac.ebi.eva.dbsnpimporter.models.Sample;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;

public class DbsnpBatchToVariantSourceProcessorTest {

    public static final int DBSNP_BUILD = 150;

    public static final int DBSNP_BATCH_ID = 12345;

    public static final String DBSNP_BATCH_NAME = "some_study";

    public static final String DBSNP_BATCH_HANDLE = "batchHandle";

    public static final String DBSNP_BATCH_HANDLE_UPP = "BATCHHANDLE";

    @Rule
    public ExpectedException exception = ExpectedException.none();

    private DbsnpBatchToVariantSourceProcessor processor;

    @Before
    public void setUp() {
        processor = new DbsnpBatchToVariantSourceProcessor(DBSNP_BUILD, DBSNP_BATCH_ID);
    }

    @Test
    public void failWithMultipleBatches() throws Exception {
        String anotherBatchName = "another_study";

        List<Sample> samples = new LinkedList<>();
        samples.add(new Sample(DBSNP_BATCH_NAME, "sample1", Sex.MALE, null, null, null));
        samples.add(new Sample(DBSNP_BATCH_NAME, "sample2", Sex.MALE, null, null, null));
        samples.add(new Sample(anotherBatchName, "sample3", Sex.MALE, null, null, null));

        DbsnpBatch dbsnpBatch = new DbsnpBatch(DBSNP_BATCH_ID, DBSNP_BATCH_HANDLE, DBSNP_BATCH_NAME, samples);

        exception.expect(IllegalArgumentException.class);
        processor.process(dbsnpBatch);
    }

    @Test
    public void failWithDuplicateSamples() throws Exception {
        Sample sample = new Sample(DBSNP_BATCH_NAME, "sample1", Sex.MALE, null, null, null);

        List<Sample> samples = new LinkedList<>();
        samples.add(sample);
        samples.add(sample);
        DbsnpBatch dbsnpBatch = new DbsnpBatch(DBSNP_BATCH_ID, DBSNP_BATCH_HANDLE, DBSNP_BATCH_NAME, samples);

        exception.expect(IllegalArgumentException.class);
        processor.process(dbsnpBatch);
    }

    @Test
    public void testProcess() throws Exception {
        Sample father = new Sample(DBSNP_BATCH_NAME, "father", Sex.MALE, null, null, null);
        Sample mother = new Sample(DBSNP_BATCH_NAME, "mother", Sex.FEMALE, null, null, null);
        Sample child = new Sample(DBSNP_BATCH_NAME, "child", Sex.MALE, father.getId(), mother.getId(), null);

        List<Sample> samples = new LinkedList<>();
        samples.add(father);
        samples.add(mother);
        samples.add(child);
        DbsnpBatch dbsnpBatch = new DbsnpBatch(DBSNP_BATCH_ID, DBSNP_BATCH_HANDLE, DBSNP_BATCH_NAME, samples);

        IVariantSource variantSource = processor.process(dbsnpBatch);

        assertEquals(DBSNP_BATCH_NAME, variantSource.getFileId());
        assertEquals(DBSNP_BATCH_HANDLE_UPP + " - " + DBSNP_BATCH_NAME, variantSource.getFileName());
        assertEquals(DBSNP_BATCH_NAME, variantSource.getStudyId());
        assertEquals(DBSNP_BATCH_HANDLE_UPP + " - " + DBSNP_BATCH_NAME, variantSource.getStudyName());

        Map<String, Integer> expectedSamplesPosition = new HashMap<>();
        expectedSamplesPosition.put(father.getName(), 0);
        expectedSamplesPosition.put(mother.getName(), 1);
        expectedSamplesPosition.put(child.getName(), 2);
        assertEquals(expectedSamplesPosition, variantSource.getSamplesPosition());

        Map<String, Object> expectedMetadata = new HashMap<>();
        expectedMetadata.put(DbsnpBatchToVariantSourceProcessor.DBSNP_BUILD_KEY, String.valueOf(DBSNP_BUILD));
        expectedMetadata.put(DbsnpBatchToVariantSourceProcessor.DBSNP_BATCH_KEY, String.valueOf(DBSNP_BATCH_ID));
        assertEquals(expectedMetadata, variantSource.getMetadata());
    }
}