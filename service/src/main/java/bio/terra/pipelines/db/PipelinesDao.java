package bio.terra.pipelines.db;

import bio.terra.pipelines.app.configuration.TspsDatabaseConfiguration;
import bio.terra.pipelines.service.model.Pipeline;
import java.util.*;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.RowMapper;
import org.springframework.stereotype.Component;

@Component
public class PipelinesDao {
  private static final RowMapper<DbPipeline> DB_PIPELINE_ROW_MAPPER =
      (rs, rowNum) ->
          new DbPipeline(
              rs.getString("pipeline_id"),
              rs.getString("display_name"),
              rs.getString("description"));

  private final JdbcTemplate tpsJdbcTemplate;

  @Autowired
  public PipelinesDao(TspsDatabaseConfiguration tspsDatabaseConfiguration) {
    this.tpsJdbcTemplate = new JdbcTemplate(tspsDatabaseConfiguration.getDataSource());
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

  private List<DbPipeline> getDbPipelines() {
    final String sql =
        """
            SELECT pipeline_id, display_name, description
            FROM pipelines
            """;

    return tpsJdbcTemplate.query(sql, DB_PIPELINE_ROW_MAPPER);
  }
}
