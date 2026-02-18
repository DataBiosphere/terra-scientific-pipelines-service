package bio.terra.pipelines.stairway.steps.datadelivery;

import bio.terra.pipelines.common.GcsFile;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.stairway.flights.datadelivery.DataDeliveryJobMapKeys;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Step to deliver pipeline output files to a GCS destination.
 *
 * <p>This step expects the following parameters in the input map: - JobMapKeys.USER_ID: the user ID
 * - DataDeliveryJobMapKeys.DESTINATION_GCS_PATH: the destination GCS path (gs://bucket/path) -
 * DataDeliveryJobMapKeys.PIPELINE_RUN_ID: the pipeline run ID to deliver outputs from
 */
public class DeliverOutputFilesToGcsStep implements Step {
  private final PipelineRunsService pipelineRunsService;
  private final PipelineInputsOutputsService pipelineInputsOutputsService;
  private final Logger logger = LoggerFactory.getLogger(DeliverOutputFilesToGcsStep.class);

  public DeliverOutputFilesToGcsStep(
      PipelineRunsService pipelineRunsService,
      PipelineInputsOutputsService pipelineInputsOutputsService) {
    this.pipelineRunsService = pipelineRunsService;
    this.pipelineInputsOutputsService = pipelineInputsOutputsService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) {
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.USER_ID,
        JobMapKeys.PIPELINE_ID,
        JobMapKeys.DOMAIN_NAME,
        DataDeliveryJobMapKeys.DESTINATION_GCS_PATH,
        DataDeliveryJobMapKeys.PIPELINE_RUN_ID);

    String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);
    GcsFile destinationGcsPath =
        inputParameters.get(DataDeliveryJobMapKeys.DESTINATION_GCS_PATH, GcsFile.class);
    UUID pipelineRunId = inputParameters.get(DataDeliveryJobMapKeys.PIPELINE_RUN_ID, UUID.class);

    logger.info(
        "Starting data delivery for pipeline run {} to destination {}",
        pipelineRunId,
        destinationGcsPath);

    PipelineRun pipelineRun = pipelineRunsService.getPipelineRun(pipelineRunId, userId);
    if (pipelineRun == null) {
      String errorMessage =
          String.format("Pipeline run %s not found for user %s", pipelineRunId, userId);
      logger.error(errorMessage);
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL, new RuntimeException(errorMessage));
    }

    try {
      pipelineInputsOutputsService.deliverOutputFilesToGcs(pipelineRun, destinationGcsPath);

      logger.info(
          "Successfully delivered output files for pipeline run {} to {}",
          pipelineRunId,
          destinationGcsPath);
      return StepResult.getStepResultSuccess();

    } catch (Exception e) {
      String errorMessage =
          String.format(
              "Failed to deliver output files for pipeline run %s to %s",
              pipelineRunId, destinationGcsPath);
      logger.error(errorMessage, e);
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_RETRY, e);
    }
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo - we'll leave the files as delivered and the delivery can be retried if
    // needed
    return StepResult.getStepResultSuccess();
  }
}
