package uk.ac.ebi.eva.dbsnpimporter.jobs.steps.processors;

import org.junit.Before;
import org.junit.Test;

import uk.ac.ebi.eva.dbsnpimporter.models.SubSnpCoreFields;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static junit.framework.TestCase.assertNotNull;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

public class MatchingAllelesFilterProcessorTest {

    private List<SubSnpCoreFields> matchingAllelesVariants;

    private List<SubSnpCoreFields> mismatchingAllelesVariants;

    private MatchingAllelesFilterProcessor filter;

    @Before
    public void setUp() throws Exception {
        filter = new MatchingAllelesFilterProcessor();
        matchingAllelesVariants = new ArrayList<>();
        matchingAllelesVariants.add(new SubSnpCoreFields(26201546, 13677177L, 1,
                                                         "NT_455866.1", 1766472L, 1766472L, 1,
                                                         "4", 91223961L, 91223961L,
                                                         "T", "T", "A", "T/A",
                                                         "NC_006091.4:g.91223961T>A", 91223961L, 91223961L, 1,
                                                         "NT_455866.1:g.1766472T>A", 1766472L, 1766472L, 1));
        matchingAllelesVariants.add(new SubSnpCoreFields(26954817, 13677177L, 1,
                                                         "NT_455866.1", 1766472L, 1766472L, 1,
                                                         "4", 91223961L, 91223961L,
                                                         "T", "T", "C", "G/A",
                                                         "NC_006091.4:g.91223961T>C", 91223961L, 91223961L, 1,
                                                         "NT_455866.1:g.1766472T>C", 1766472L, 1766472L, 1));
        matchingAllelesVariants.add(new SubSnpCoreFields(26963037, 13677177L, 1,
                                                         "NT_455866.1", 1766472L, 1766472L, 1,
                                                         "4", 91223961L, 91223961L,
                                                         "T", "T", "A", "T/A",
                                                         "NC_006091.4:g.91223961T>A", 91223961L, 91223961L, 1,
                                                         "NT_455866.1:g.1766472T>A", 1766472L, 1766472L, 1));

        mismatchingAllelesVariants = new ArrayList<>();
        mismatchingAllelesVariants.add(new SubSnpCoreFields(26201546, 13677177L, 1,
                                                            "NT_455866.1", 1766472L, 1766472L, 1,
                                                            "4", 91223961L, 91223961L,
                                                            "T", "T", "C", "T/A",
                                                            "NC_006091.4:g.91223961T>C", 91223961L, 91223961L, 1,
                                                            "NT_455866.1:g.1766472T>C", 1766472L, 1766472L, 1));
        mismatchingAllelesVariants.add(new SubSnpCoreFields(26954817, 13677177L, 1,
                                                            "NT_455866.1", 1766472L, 1766472L, 1,
                                                            "4", 91223961L, 91223961L,
                                                            "T", "T", "A", "G/A",
                                                            "NC_006091.4:g.91223961T>A", 91223961L, 91223961L, 1,
                                                            "NT_455866.1:g.1766472T>A", 1766472L, 1766472L, 1));
        mismatchingAllelesVariants.add(new SubSnpCoreFields(26963037, 13677177L, 1,
                                                            "NT_455866.1", 1766472L, 1766472L, 1,
                                                            "4", 91223961L, 91223961L,
                                                            "T", "T", "C", "T/A",
                                                            "NC_006091.4:g.91223961T>C", 91223961L, 91223961L, 1,
                                                            "NT_455866.1:g.1766472T>C", 1766472L, 1766472L, 1));
    }

    @Test
    public void removeMismatchingVariant() throws Exception {
        assertNull(filter.process(mismatchingAllelesVariants.get(0)));
    }

    @Test
    public void keepMatchingVariant() throws Exception {
        assertNotNull(filter.process(matchingAllelesVariants.get(0)));
    }

    @Test
    public void mismatchingVariantsMustBeRemoved() throws Exception {
        List<SubSnpCoreFields> mixedVariants = new ArrayList<>();
        mixedVariants.addAll(mismatchingAllelesVariants);
        mixedVariants.addAll(matchingAllelesVariants);
        Collections.shuffle(mixedVariants);

        Set<SubSnpCoreFields> filteredVariants = new HashSet<>();
        for (SubSnpCoreFields subSnp : mixedVariants) {
            if (filter.process(subSnp) != null) {
                filteredVariants.add(subSnp);
            }
        }
        assertEquals(new HashSet<>(matchingAllelesVariants), filteredVariants);
    }
}