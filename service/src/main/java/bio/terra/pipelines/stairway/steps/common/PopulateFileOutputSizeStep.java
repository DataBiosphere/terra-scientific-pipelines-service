package bio.terra.pipelines.stairway.steps.common;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.model.Pipeline;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import io.sentry.Sentry;
import io.sentry.SentryLevel;
import java.util.Map;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Step to populate the file sizes of pipeline outputs in the working map.
 *
 * <p>This step expects the JobMapKeys.PIPELINE_KEY in the input parameters and
 * ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS in the working map. It will write the output file sizes
 * to ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE in the working map.
 *
 * <p>If there is an error populating the file sizes, this step will log the error and continue
 * without populating the file sizes, since we don't want to fail the entire pipeline run if we
 * can't get the file sizes.
 */
public class PopulateFileOutputSizeStep implements Step {
  private final PipelinesService pipelinesService;
  private final PipelineInputsOutputsService pipelineInputsOutputsService;
  private final Logger logger = LoggerFactory.getLogger(PopulateFileOutputSizeStep.class);

  public PopulateFileOutputSizeStep(
      PipelinesService pipelinesService,
      PipelineInputsOutputsService pipelineInputsOutputsService) {
    this.pipelinesService = pipelinesService;
    this.pipelineInputsOutputsService = pipelineInputsOutputsService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) {
    try {
      // validate and extract parameters from input map
      var inputParameters = flightContext.getInputParameters();
      FlightUtils.validateRequiredEntries(inputParameters, JobMapKeys.PIPELINE_KEY);
      String pipelineKey = inputParameters.get(JobMapKeys.PIPELINE_KEY, String.class);

      // validate and extract parameters from working map
      var workingMap = flightContext.getWorkingMap();
      FlightUtils.validateRequiredEntries(workingMap, ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS);
      Map<String, String> outputsMap =
          workingMap.get(ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS, Map.class);

      Pipeline pipeline = pipelinesService.getPipelineByKey(pipelineKey);
      Map<String, Long> outputFileSizes =
          pipelineInputsOutputsService.getPipelineOutputsFileSize(pipeline, outputsMap);

      logger.info("Retrieved file sizes for pipeline outputs");
      workingMap.put(ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE, outputFileSizes);
    } catch (Exception e) {
      // log the error and continue without populating the file sizes as we don't want to fail
      // the entire pipeline run if we can't get the file sizes
      logger.error(
          "Error populating file output sizes for pipeline run, continuing without them. The step will be marked as success.",
          e);
      Sentry.captureMessage(
          "Error populating file output sizes for flight id %s. Error details : %s"
              .formatted(flightContext.getFlightId(), e.getMessage()),
          SentryLevel.ERROR);
    }

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo
    return StepResult.getStepResultSuccess();
  }
}
