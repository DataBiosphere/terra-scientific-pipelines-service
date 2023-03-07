package bio.terra.pipelines.db;

import bio.terra.pipelines.app.configuration.TspsDatabaseConfiguration;
import bio.terra.pipelines.service.model.Pipeline;
import java.util.*;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.jdbc.core.RowMapper;
import org.springframework.jdbc.core.namedparam.MapSqlParameterSource;
import org.springframework.jdbc.core.namedparam.NamedParameterJdbcTemplate;
import org.springframework.stereotype.Component;

@Component
public class PipelinesDao {
  private static final RowMapper<DbPipeline> DB_PIPELINE_ROW_MAPPER =
      (rs, rowNum) ->
          new DbPipeline(
              rs.getString("pipeline_id"),
              rs.getString("display_name"),
              rs.getString("description"));

  private final NamedParameterJdbcTemplate jdbcTemplate;

  @Autowired
  public PipelinesDao(TspsDatabaseConfiguration tspsDatabaseConfiguration) {
    this.jdbcTemplate = new NamedParameterJdbcTemplate(tspsDatabaseConfiguration.getDataSource());
  }

  /**
   * Check that a pipeline exists in the database exactly once
   *
   * @param pipelineId
   * @return boolean
   */
  public boolean checkPipelineExists(String pipelineId) {
    Integer dbPipelineCount = countDbPipeline(pipelineId);

    return dbPipelineCount == 1;
  }

  /**
   * Return all available Pipelines
   *
   * @return List of Pipeline objects
   */
  public List<Pipeline> getPipelines() {
    List<Pipeline> pipelineList = new ArrayList<>();

    List<DbPipeline> dbPipelineList = getDbPipelines();

    for (DbPipeline dbPipeline : dbPipelineList) {
      pipelineList.add(Pipeline.fromDb(dbPipeline));
    }

    return pipelineList;
  }

  private Integer countDbPipeline(String pipelineId) {
    final String sql =
        """
                    SELECT COUNT(pipeline_id)
                    FROM pipelines
                    WHERE pipeline_id = :pipelineId
                    """;
    MapSqlParameterSource params = new MapSqlParameterSource().addValue("pipelineId", pipelineId);

    return jdbcTemplate.queryForObject(sql, params, Integer.class);
  }

  private List<DbPipeline> getDbPipelines() {
    final String sql =
        """
            SELECT pipeline_id, display_name, description
            FROM pipelines
            """;

    return jdbcTemplate.query(sql, DB_PIPELINE_ROW_MAPPER);
  }
}
