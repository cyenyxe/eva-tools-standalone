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
package uk.ac.ebi.eva.dbsnpimporter.io;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;

import uk.ac.ebi.eva.dbsnpimporter.models.PopulationFrequencies;

import java.io.IOException;
import java.util.List;

public class FrequencyInfoParser {

    private final ObjectMapper objectMapper;

    private final TypeReference<List<PopulationFrequencies>> typeReference;

    public FrequencyInfoParser() {
        objectMapper = new ObjectMapper();
        typeReference = new TypeReference<List<PopulationFrequencies>>() { };
    }

    public List<PopulationFrequencies> parse(String frequencyInfo) throws IOException {
        return objectMapper.readValue(frequencyInfo, typeReference);
    }
}