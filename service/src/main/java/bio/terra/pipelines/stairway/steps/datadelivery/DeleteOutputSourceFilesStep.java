package bio.terra.pipelines.stairway.steps.datadelivery;

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
import io.sentry.Sentry;
import io.sentry.SentryLevel;
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
public class DeleteOutputSourceFilesStep implements Step {
  private final PipelineRunsService pipelineRunsService;
  private final PipelineInputsOutputsService pipelineInputsOutputsService;
  private final Logger logger = LoggerFactory.getLogger(DeliverOutputFilesToGcsStep.class);

  public DeleteOutputSourceFilesStep(
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
        JobMapKeys.DOMAIN_NAME,
        DataDeliveryJobMapKeys.PIPELINE_RUN_ID);

    String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);
    UUID pipelineRunId = inputParameters.get(DataDeliveryJobMapKeys.PIPELINE_RUN_ID, UUID.class);

    logger.info(
        "Starting deletion of source output files for pipeline run {}, after successful delivery to destination",
        pipelineRunId);

    PipelineRun pipelineRun = pipelineRunsService.getPipelineRun(pipelineRunId, userId);
    if (pipelineRun == null) {
      String errorMessage =
          String.format("Pipeline run %s not found for user %s", pipelineRunId, userId);
      logger.error(errorMessage);
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL, new RuntimeException(errorMessage));
    }

    try {
      pipelineInputsOutputsService.deleteOutputSourcesFiles(pipelineRun);

      logger.info(
          "Successfully deleted source output files for pipeline run {} after successful delivery to destination",
          pipelineRunId);
      return StepResult.getStepResultSuccess();

    } catch (Exception e) {
      String errorMessage =
          String.format("Failed to delete output source files for pipeline run %s", pipelineRunId);
      logger.error(errorMessage, e);
      Sentry.captureMessage(errorMessage, SentryLevel.ERROR);
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
