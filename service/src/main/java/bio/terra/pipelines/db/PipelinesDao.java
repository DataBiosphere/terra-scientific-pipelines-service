package bio.terra.pipelines.db;

import bio.terra.pipelines.app.configuration.TspsDatabaseConfiguration;
import bio.terra.pipelines.service.model.Pipeline;
import java.util.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.RowMapper;
import org.springframework.stereotype.Component;

@Component
public class PipelinesDao {
  private static final RowMapper<DbPipeline> DB_PIPELINE_ROW_MAPPER =
      (rs, rowNum) -> {
        return new DbPipeline(
            rs.getString("pipeline_id"), rs.getString("display_name"), rs.getString("description"));
      };

  private final Logger logger = LoggerFactory.getLogger(PipelinesDao.class);
  private final JdbcTemplate tpsJdbcTemplate;

  @Autowired
  public PipelinesDao(TspsDatabaseConfiguration tspsDatabaseConfiguration) {
    this.tpsJdbcTemplate = new JdbcTemplate(tspsDatabaseConfiguration.getDataSource());
  }

  //  @ReadTransaction
  //  public Pipeline getPao(UUID objectId) {
  //    DbPipeline dbPao = getDbPao(objectId);
  //    Map<String, PolicyInputs> attributeSetMap =
  //        getAttributeSets(List.of(dbPao.attributeSetId(), dbPao.effectiveSetId()));
  //    return Pipeline.fromDb(dbPao, attributeSetMap);
  //  }

  // -- Graph Walk Methods --
  // These methods are intentionally without transaction annotations. They are used by the policy
  // update process. That process may do multiple reads of the database followed by a big update.
  // The transaction is managed outside of the DAO.

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
