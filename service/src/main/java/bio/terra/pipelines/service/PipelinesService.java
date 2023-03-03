package bio.terra.pipelines.service;

import bio.terra.common.exception.NotFoundException;
import bio.terra.pipelines.db.PipelinesDao;
import bio.terra.pipelines.service.model.Pipeline;
import java.util.List;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

/** The Pipelines Service manages information about the service's available Scientific Pipelines. */
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

  public boolean validatePipeline(String pipelineId) {
    logger.info("Validate pipeline");

    try {
      pipelinesDao.getPipeline(pipelineId);
      return true;
    } catch (NotFoundException e) {
      return false;
    }
  }
}
