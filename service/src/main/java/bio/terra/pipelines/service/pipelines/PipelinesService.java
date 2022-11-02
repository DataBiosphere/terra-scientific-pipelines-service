package bio.terra.pipelines.service.pipelines;

import bio.terra.pipelines.db.PipelinesDao;
import bio.terra.pipelines.service.pipelines.model.Pipeline;
import java.util.List;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

/**
 * The Pipelines Service does all processing of the Policy Attribute Objects. Those are the objects
 * that represent policies on objects in other components of Terra. It manages the graph (DAG) of
 * objects that depend on the policies of other objects.
 */
@Component
public class PipelinesService {
  private static final Logger logger = LoggerFactory.getLogger(PipelinesService.class);

  private final PipelinesDao pipelinesDao;

  @Autowired
  public PipelinesService(PipelinesDao pipelinesDao) {
    this.pipelinesDao = pipelinesDao;
  }

  public List<Pipeline> getPipelines() {
    logger.info("Get all Pipelines");
    return pipelinesDao.getPipelines();
  }
}
