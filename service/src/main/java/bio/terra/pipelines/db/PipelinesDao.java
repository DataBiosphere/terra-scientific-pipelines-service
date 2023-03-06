package bio.terra.pipelines.db;

import bio.terra.common.exception.NotFoundException;
import bio.terra.pipelines.app.configuration.TspsDatabaseConfiguration;
import bio.terra.pipelines.service.model.Pipeline;
import java.util.*;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.EmptyResultDataAccessException;
import org.springframework.dao.support.DataAccessUtils;
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
   * Return the matching pipeline, given the pipelineId and version
   *
   * @param pipelineId
   * @return Pipeline object
   */
  public Pipeline getPipeline(String pipelineId) {
    DbPipeline dbPipeline =
        getDbPipelineIfExists(pipelineId)
            .orElseThrow(
                () -> new NotFoundException(String.format("Pipeline %s not found.", pipelineId)));

    return Pipeline.fromDb(dbPipeline);
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

  private Optional<DbPipeline> getDbPipelineIfExists(String pipelineId) {
    final String sql =
        """
                    SELECT pipeline_id, display_name, description
                    FROM pipelines
                    WHERE pipeline_id = :pipelineId
                    """;
    MapSqlParameterSource params = new MapSqlParameterSource().addValue("pipelineId", pipelineId);

    try {
      DbPipeline result =
          DataAccessUtils.requiredSingleResult(
              jdbcTemplate.query(sql, params, DB_PIPELINE_ROW_MAPPER));
      return Optional.of(result);
    } catch (EmptyResultDataAccessException e) {
      return Optional.empty();
    }
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
