package bio.terra.pipelines.stairway.steps.datadelivery;

import bio.terra.pipelines.common.GcsFile;
import bio.terra.pipelines.common.utils.DataDeliveryStatusEnum;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.DataDeliveryService;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.stairway.flights.datadelivery.DataDeliveryJobMapKeys;
import bio.terra.pipelines.stairway.steps.exception.InternalStepException;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Step to create a data delivery record in the database with status RUNNING.
 *
 * <p>This step expects the following parameters in the input map: - JobMapKeys.USER_ID: the user ID
 * - DataDeliveryJobMapKeys.DESTINATION_GCS_PATH: the destination GCS path (gs://bucket/path) -
 * DataDeliveryJobMapKeys.PIPELINE_RUN_ID: the pipeline run ID to deliver outputs from
 */
public class CreateDataDeliveryRecordStep implements Step {
  private final PipelineRunsService pipelineRunsService;
  private final DataDeliveryService dataDeliveryService;
  private final Logger logger = LoggerFactory.getLogger(CreateDataDeliveryRecordStep.class);

  public CreateDataDeliveryRecordStep(
      PipelineRunsService pipelineRunsService, DataDeliveryService dataDeliveryService) {
    this.pipelineRunsService = pipelineRunsService;
    this.dataDeliveryService = dataDeliveryService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) {
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.USER_ID,
        DataDeliveryJobMapKeys.DESTINATION_GCS_PATH,
        DataDeliveryJobMapKeys.PIPELINE_RUN_ID);

    String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);
    GcsFile destinationGcsPath =
        inputParameters.get(DataDeliveryJobMapKeys.DESTINATION_GCS_PATH, GcsFile.class);
    UUID pipelineRunId = inputParameters.get(DataDeliveryJobMapKeys.PIPELINE_RUN_ID, UUID.class);

    PipelineRun pipelineRun = pipelineRunsService.getPipelineRun(pipelineRunId, userId);
    if (pipelineRun == null) {
      String errorMessage =
          String.format("Pipeline run %s not found for user %s", pipelineRunId, userId);
      logger.error(errorMessage);
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL, new InternalStepException("Data delivery failed."));
    }

    logger.info(
        "Creating data delivery record for pipeline run {} to destination {}",
        pipelineRunId,
        destinationGcsPath);

    dataDeliveryService.createDataDelivery(
        pipelineRun.getId(),
        UUID.fromString(flightContext.getFlightId()),
        DataDeliveryStatusEnum.RUNNING,
        destinationGcsPath);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    var inputParameters = flightContext.getInputParameters();
    String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);
    UUID pipelineRunId = inputParameters.get(DataDeliveryJobMapKeys.PIPELINE_RUN_ID, UUID.class);

    PipelineRun pipelineRun = pipelineRunsService.getPipelineRun(pipelineRunId, userId);
    if (pipelineRun != null) {
      dataDeliveryService.updateDataDeliveryStatus(
          pipelineRun.getId(), DataDeliveryStatusEnum.FAILED);
    }
    return StepResult.getStepResultSuccess();
  }
}
