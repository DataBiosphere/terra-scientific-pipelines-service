package bio.terra.pipelines.stairway.steps.datadelivery;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.stairway.flights.datadelivery.DataDeliveryJobMapKeys;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Step to deliver pipeline output files to a GCS destination.
 *
 * <p>This step expects the following parameters in the input map: - JobMapKeys.USER_ID: the user ID
 * - DataDeliveryJobMapKeys.DESTINATION_GCS_PATH: the destination GCS path (gs://bucket/path) -
 * DataDeliveryJobMapKeys.PIPELINE_RUN_ID: the pipeline run ID to get outputs from
 */
public class DeleteOutputSourceFilesStep implements Step {
  private final PipelineRunsService pipelineRunsService;
  private final PipelineInputsOutputsService pipelineInputsOutputsService;
  private final PipelinesService pipelinesService;
  private final Logger logger = LoggerFactory.getLogger(DeliverOutputFilesToGcsStep.class);

  public DeleteOutputSourceFilesStep(
      PipelineRunsService pipelineRunsService,
      PipelineInputsOutputsService pipelineInputsOutputsService,
      PipelinesService pipelinesService) {
    this.pipelineRunsService = pipelineRunsService;
    this.pipelineInputsOutputsService = pipelineInputsOutputsService;
    this.pipelinesService = pipelinesService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) {
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.USER_ID,
        JobMapKeys.DOMAIN_NAME,
        DataDeliveryJobMapKeys.DESTINATION_GCS_PATH,
        DataDeliveryJobMapKeys.PIPELINE_RUN_ID);

    UUID pipelineRunId = inputParameters.get(DataDeliveryJobMapKeys.PIPELINE_RUN_ID, UUID.class);
    String destinationGcsPath =
        inputParameters.get(DataDeliveryJobMapKeys.DESTINATION_GCS_PATH, String.class);

    logger.info(
        "(no-op for now) Successfully deleted original output files for pipeline run {} to {}",
        pipelineRunId,
        destinationGcsPath);
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    logger.info(
        "Undo step for public class DeleteOutputSourceFilesStep implements Step {\n - no-op for now");
    return StepResult.getStepResultSuccess();
  }
}
