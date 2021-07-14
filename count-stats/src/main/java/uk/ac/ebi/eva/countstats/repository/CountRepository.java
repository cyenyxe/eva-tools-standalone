package uk.ac.ebi.eva.countstats.repository;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.jdbc.core.namedparam.MapSqlParameterSource;
import org.springframework.jdbc.core.namedparam.NamedParameterJdbcOperations;
import org.springframework.jdbc.core.namedparam.SqlParameterSource;
import org.springframework.stereotype.Repository;
import uk.ac.ebi.eva.countstats.model.CountDto;

import java.math.BigInteger;
import java.sql.Types;
import java.util.List;

@Repository
public class CountRepository {
    @Autowired
    private NamedParameterJdbcOperations namedParameterJdbcOperations;

    public int saveCount(CountDto countDto) {
        String sql = "insert into process_count_metric values (:process, :identifier, :metric, :count)";

        SqlParameterSource sqlParameterSource = new MapSqlParameterSource("process", countDto.getProcess())
                .addValue("identifier", countDto.getIdentifier(), Types.OTHER)
                .addValue("metric", countDto.getMetric())
                .addValue("count", countDto.getCount(), Types.BIGINT);

        return namedParameterJdbcOperations.update(sql, sqlParameterSource);
    }

    public List<CountDto> getAllCounts() {
        String sql = "select * from process_count_metric";

        return namedParameterJdbcOperations.query(sql, (rs, rownu) -> new CountDto(rs.getString("process"),
                rs.getString("identifier"), rs.getString("metric"), rs.getObject("count", BigInteger.class))
        );
    }
}
