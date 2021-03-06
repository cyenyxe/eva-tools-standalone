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
package uk.ac.ebi.eva.dbsnpimporter.io.readers;

import org.springframework.jdbc.core.RowMapper;

import uk.ac.ebi.eva.dbsnpimporter.models.LocusType;
import uk.ac.ebi.eva.dbsnpimporter.models.Orientation;
import uk.ac.ebi.eva.dbsnpimporter.models.SubSnpCoreFields;

import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Maps the database fields that correspond to an SS ID along with its contig and (optionally) chromosome coordinates.
 */
public class SubSnpCoreFieldsRowMapper implements RowMapper<SubSnpCoreFields> {

    public static final String SUBSNP_ID_COLUMN = "ss_id";

    public static final String REFSNP_ID_COLUMN = "rs_id";

    public static final String LOC_TYPE_COLUMN = "loc_type";

    public static final String CONTIG_NAME_COLUMN = "contig_name";

    public static final String CONTIG_START_COLUMN = "contig_start";

    public static final String CONTIG_END_COLUMN = "contig_end";

    public static final String CHROMOSOME_COLUMN = "chromosome";

    public static final String CHROMOSOME_START_COLUMN = "chromosome_start";

    public static final String CHROMOSOME_END_COLUMN = "chromosome_end";

    public static final String SNP_ORIENTATION_COLUMN = "snp_orientation";

    public static final String CONTIG_ORIENTATION_COLUMN = "contig_orientation";

    public static final String HGVS_C_STRING = "hgvs_c_string";

    public static final String HGVS_C_START = "hgvs_c_start";

    public static final String HGVS_C_STOP = "hgvs_c_stop";

    public static final String REFERENCE_C = "reference_c";

    public static final String HGVS_C_ORIENTATION = "hgvs_c_orientation";

    public static final String HGVS_T_STRING = "hgvs_t_string";

    public static final String HGVS_T_START = "hgvs_t_start";

    public static final String HGVS_T_STOP = "hgvs_t_stop";

    public static final String REFERENCE_T = "reference_t";

    public static final String HGVS_T_ORIENTATION = "hgvs_t_orientation";

    public static final String ALTERNATE = "alternate";

    public static final String ALLELES = "alleles";

    public static final String SUBSNP_ORIENTATION_COLUMN = "subsnp_orientation";

    public static final String BATCH_COLUMN = "batch_name";

    public static final String LOAD_ORDER_COLUMN = "load_order";

    public static final String GENOTYPES_COLUMN = "genotypes_string";

    public static final String FREQUENCIES_COLUMN = "freq_info";

    private ResultSet resultSet;

    /**
     * Maps ResultSet to SubSnpCoreFields.
     *
     * It makes getObject (instead of getInt or getString) for those that are nullable.
     *
     * The casts are safe because the DB types are integers. The types Long and BigDecimal are introduced by the query
     * and it won't change the values more that +-1.
     */
    @Override
    public SubSnpCoreFields mapRow(ResultSet resultSet, int rowNum) throws SQLException {
        this.resultSet = resultSet;
        return new SubSnpCoreFields(
                resultSet.getLong(SUBSNP_ID_COLUMN),
                Orientation.getOrientation(resultSet.getInt(SUBSNP_ORIENTATION_COLUMN)),
                getAsLong(REFSNP_ID_COLUMN),
                Orientation.getOrientation(resultSet.getInt(SNP_ORIENTATION_COLUMN)),
                resultSet.getString(CONTIG_NAME_COLUMN),
                resultSet.getLong(CONTIG_START_COLUMN),
                resultSet.getLong(CONTIG_END_COLUMN),
                Orientation.getOrientation(resultSet.getInt(CONTIG_ORIENTATION_COLUMN)),
                LocusType.getLocusType(resultSet.getInt(LOC_TYPE_COLUMN)),
                resultSet.getString(CHROMOSOME_COLUMN),
                getAsLong(CHROMOSOME_START_COLUMN),
                getAsLong(CHROMOSOME_END_COLUMN),
                resultSet.getString(REFERENCE_C),
                resultSet.getString(REFERENCE_T),
                resultSet.getString(ALTERNATE),
                resultSet.getString(ALLELES),
                resultSet.getString(HGVS_C_STRING),
                getAsLong(HGVS_C_START),
                getAsLong(HGVS_C_STOP),
                Orientation.getOrientation(resultSet.getInt(HGVS_C_ORIENTATION)),
                resultSet.getString(HGVS_T_STRING),
                getAsLong(HGVS_T_START),
                getAsLong(HGVS_T_STOP),
                Orientation.getOrientation(resultSet.getInt(HGVS_T_ORIENTATION)),
                resultSet.getString(GENOTYPES_COLUMN),
                resultSet.getString(FREQUENCIES_COLUMN),
                resultSet.getString(BATCH_COLUMN));
    }

    private Long getAsLong(String column) throws SQLException {
        Object object = resultSet.getObject(column);
        if (object == null) {
            return null;
        } else {
            if (object instanceof Number) {
                return ((Number) object).longValue();
            } else {
                int columnIndex = resultSet.findColumn(column);
                int columnType = resultSet.getMetaData().getColumnType(columnIndex);
                throw new IllegalArgumentException("Can not convert column '" + column + "' of type " + columnType
                                                           + " (see java.sql.Types) to Number.");
            }
        }
    }
}
