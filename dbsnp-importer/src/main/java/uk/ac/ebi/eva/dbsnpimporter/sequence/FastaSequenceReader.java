/*
 * Copyright 2017 EMBL - European Bioinformatics Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package uk.ac.ebi.eva.dbsnpimporter.sequence;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

import java.nio.file.Path;

/**
 * Implementation of SequenceReader for indexed fasta files
 */
public class FastaSequenceReader implements SequenceReader {

    private final ReferenceSequenceFile fastaSequenceFile;

    public FastaSequenceReader(Path fastaFile) {
        fastaSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(fastaFile, true);
    }

    @Override
    public String getSequence(String contig, long start, long end) throws ReadSequenceException {
        if (end < start) {
            throw new ReadSequenceException("'end' should be greater or equal than 'start'");
        } else if (start < 1) {
            throw new ReadSequenceException("'start' and 'end' should be positive integers");
        }

        try {
            return fastaSequenceFile.getSubsequenceAt(contig, start, end).getBaseString();
        } catch (SAMException e) {
            throw new ReadSequenceException(e);
        }
    }

    @Override
    public void close() throws Exception {
        fastaSequenceFile.close();
    }
}