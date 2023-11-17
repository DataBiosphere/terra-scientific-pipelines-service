package bio.terra.pipelines.service;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.dependencies.stairway.StairwayJobBuilder;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.stairway.GetPipelineFlight;
import java.util.List;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

/** The Pipelines Service manages information about the service's available Scientific Pipelines. */
@Component
public class PipelinesService {
  private static final Logger logger = LoggerFactory.getLogger(PipelinesService.class);

  private final PipelinesRepository pipelinesRepository;
  private final StairwayJobService stairwayJobService;

  @Autowired
  public PipelinesService(
      PipelinesRepository pipelinesRepository, StairwayJobService stairwayJobService) {
    this.pipelinesRepository = pipelinesRepository;
    this.stairwayJobService = stairwayJobService;
  }

  public List<Pipeline> getPipelines() {
    logger.info("Get all Pipelines");
    return pipelinesRepository.findAll();
  }

  public Pipeline getImputationPipelineViaFlight() {
    logger.info("Get the imputation Pipeline via flight - a toy flight example");
    StairwayJobBuilder stairwayJobBuilder =
        stairwayJobService.newJob().flightClass(GetPipelineFlight.class);
    return stairwayJobBuilder.submitAndWait(Pipeline.class);
  }

  public Pipeline getPipeline(String pipelineId) {
    logger.info("Get a specific pipeline for pipelineId {}".format(pipelineId));
    return pipelinesRepository.findByPipelineId(pipelineId);
  }

  public boolean pipelineExists(String pipelineId) {
    return pipelinesRepository.existsByPipelineId(pipelineId);
  }
}
