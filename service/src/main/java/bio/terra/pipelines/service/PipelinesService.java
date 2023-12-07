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

  public Pipeline getImputationPipelineViaFlight(String pipelineId) {
    StairwayJobBuilder stairwayJobBuilder =
        stairwayJobService.newJob().flightClass(GetPipelineFlight.class).pipelineId(pipelineId);

    // normally we'd submit the job and use a different endpoint to retrieve the result, but since
    // this is a quick toy example, we submit it, wait a second for it to complete, and then
    // retrieve the result
    String jobId = stairwayJobBuilder.submit();

    // wait for job to complete
    try {
      Thread.sleep(1000);
    } catch (InterruptedException e) {
      logger.error("Interrupted while waiting for job to complete", e);
      Thread.currentThread().interrupt();
    }

    logger.info("Retrieving flight result (jobId: {})", jobId);
    return stairwayJobService.retrieveJobResult(jobId, Pipeline.class).getResult();
  }

  public Pipeline getPipeline(String pipelineId) {
    logger.info("Get a specific pipeline for pipelineId {}", pipelineId);
    Pipeline dbResult = pipelinesRepository.findByPipelineId(pipelineId);
    if (dbResult == null) {
      throw new IllegalArgumentException(
          String.format("Pipeline not found for pipelineId %s", pipelineId));
    }
    return dbResult;
  }

  public boolean pipelineExists(String pipelineId) {
    return pipelinesRepository.existsByPipelineId(pipelineId);
  }
}
