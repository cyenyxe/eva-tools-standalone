/*
 *
 *  * Copyright 2016 EMBL - European Bioinformatics Institute
 *  *
 *  * Licensed under the Apache License, Version 2.0 (the "License");
 *  * you may not use this file except in compliance with the License.
 *  * You may obtain a copy of the License at
 *  *
 *  *      http://www.apache.org/licenses/LICENSE-2.0
 *  *
 *  * Unless required by applicable law or agreed to in writing, software
 *  * distributed under the License is distributed on an "AS IS" BASIS,
 *  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  * See the License for the specific language governing permissions and
 *  * limitations under the License.
 *
 */

package uk.ac.ebi.eva.vcfdump.regionutils;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.eva.commons.core.models.Region;
import uk.ac.ebi.eva.commons.mongodb.services.VariantWithSamplesAndAnnotationsService;
import uk.ac.ebi.eva.vcfdump.QueryParams;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public class RegionFactory {

    private static final Logger logger = LoggerFactory.getLogger(RegionFactory.class);

    private int windowSize;

    private final VariantWithSamplesAndAnnotationsService variantService;

    public RegionFactory(int windowSize, VariantWithSamplesAndAnnotationsService variantService) {
        this.windowSize = windowSize;
        this.variantService = variantService;
    }

    public List<Region> getRegionsForChromosome(String chromosome, QueryParams query) {
        String regionFilter = query.getRegion();
        if (regionFilter == null || regionFilter.isEmpty() || isChromosomeInRegionFilterWithNoCoordinates(chromosome,
                                                                                                          regionFilter)) {
            // if there are no region filter or no chromosome coordinates in the filter, we need to get the min and max variant start from mongo
            Long minStart = variantService.findChromosomeLowestReportedCoordinate(chromosome, query.getStudies());
            if (minStart == null) {
                return Collections.EMPTY_LIST;
            } else {
                long maxStart = variantService.findChromosomeHighestReportedCoordinate(chromosome, query.getStudies());;
                logger.debug("Chromosome {} maxStart: {}", chromosome, maxStart);
                logger.debug("Chromosome {} minStart: {}", chromosome, minStart);
                return divideChromosomeInChunks(chromosome, minStart, maxStart);
            }
        } else {
            List<Region> chromosomeRegionsFromQuery =
                    Region.parseRegions(regionFilter).stream().filter(r -> r.getChromosome().equals(chromosome))
                                                     .collect(new IntersectingRegionsMerger());

            String commaSeparatedRegionList = chromosomeRegionsFromQuery.stream().map(Region::toString)
                                                                        .collect(Collectors.joining(", "));
            logger.debug("Chromosome {} regions from query: {}", chromosome, commaSeparatedRegionList);

            return divideRegionListInChunks(chromosomeRegionsFromQuery);

        }
    }

    private boolean isChromosomeInRegionFilterWithNoCoordinates(String chromosome, String regionFilter) {
        return Arrays.asList(regionFilter.split(",")).stream()
                     .anyMatch(regionString -> regionString.equals(chromosome));
    }

    public List<Region> divideChromosomeInChunks(String chromosome, long chromosomeMinStart, long chromosomeMaxStart) {
        List<Region> regions = divideRegionInChunks(chromosome, chromosomeMinStart, chromosomeMaxStart);
        logger.debug("Number of regions in chromosome{}: {}", chromosome, regions.size());
        if (!regions.isEmpty()) {
            logger.debug("First region: {}", regions.get(0));
            logger.debug("Last region: {}", regions.get(regions.size() - 1));
        }

        return regions;
    }


    private List<Region> divideRegionListInChunks(List<Region> regionsFromQuery) {
        List<Region> regions = new ArrayList<>();
        for (Region region : regionsFromQuery) {
            regions.addAll(divideRegionInChunks(region.getChromosome(), region.getStart(), region.getEnd()));
        }
        return regions;
    }

    public List<Region> divideRegionInChunks(String chromosome, long minStart, long maxStart) {
        // we are using long instead of int to avoid overflowing if maxStart is MAX_INT. The casting to int will always work, because the
        // maximum value will be MAX_INT
        if (minStart == -1) {
            return Collections.EMPTY_LIST;
        } else {
            List<Region> regions = new ArrayList<>();
            long nextStart = minStart;
            while (nextStart <= maxStart) {
                long end = Math.min(nextStart + windowSize, maxStart + 1);
                regions.add(new Region(chromosome, nextStart, (end - 1)));
                nextStart = end;
            }
            return regions;
        }
    }

}
