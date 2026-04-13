package bio.terra.pipelines.stairway.steps.datadelivery;

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
 * Step to update the data delivery record status to SUCCEEDED after a successful delivery.
 *
 * <p>This step expects the following parameters in the input map: - JobMapKeys.USER_ID: the user ID
 * - DataDeliveryJobMapKeys.PIPELINE_RUN_ID: the pipeline run ID whose delivery record to update
 */
public class UpdateDataDeliveryStatusStep implements Step {
  private final PipelineRunsService pipelineRunsService;
  private final DataDeliveryService dataDeliveryService;
  private final Logger logger = LoggerFactory.getLogger(UpdateDataDeliveryStatusStep.class);

  public UpdateDataDeliveryStatusStep(
      PipelineRunsService pipelineRunsService, DataDeliveryService dataDeliveryService) {
    this.pipelineRunsService = pipelineRunsService;
    this.dataDeliveryService = dataDeliveryService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) {
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters, JobMapKeys.USER_ID, DataDeliveryJobMapKeys.PIPELINE_RUN_ID);

    String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);
    UUID pipelineRunId = inputParameters.get(DataDeliveryJobMapKeys.PIPELINE_RUN_ID, UUID.class);

    PipelineRun pipelineRun = pipelineRunsService.getPipelineRun(pipelineRunId, userId);
    if (pipelineRun == null) {
      String errorMessage =
          String.format("Pipeline run %s not found for user %s", pipelineRunId, userId);
      logger.error(errorMessage);
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL, new InternalStepException("Data delivery failed."));
    }

    logger.info("Marking data delivery record as SUCCEEDED for pipeline run {}", pipelineRunId);

    dataDeliveryService.updateDataDeliveryStatus(
        pipelineRun.getId(), DataDeliveryStatusEnum.SUCCEEDED);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo - the files have already been successfully delivered
    return StepResult.getStepResultSuccess();
  }
}
