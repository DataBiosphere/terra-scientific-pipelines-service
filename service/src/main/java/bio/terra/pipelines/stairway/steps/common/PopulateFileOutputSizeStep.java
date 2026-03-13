package bio.terra.pipelines.stairway.steps.common;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.stairway.flights.imputation.ImputationJobMapKeys;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import java.util.Map;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

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
    // validate and extract parameters from input map
    var inputParameters = flightContext.getInputParameters();
    FlightUtils.validateRequiredEntries(inputParameters, JobMapKeys.PIPELINE_ID);
    Long pipelineId = inputParameters.get(JobMapKeys.PIPELINE_ID, Long.class);

    // validate and extract parameters from working map
    var workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(workingMap, ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS);
    Map<String, String> outputsMap =
        workingMap.get(ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS, Map.class);

    Pipeline pipeline = pipelinesService.getPipelineById(pipelineId);
    Map<String, Long> outputFileSizes =
        pipelineInputsOutputsService.getPipelineOutputsFileSize(pipeline, outputsMap);

    logger.info("Retrieved file sizes for pipeline outputs");
    workingMap.put(ImputationJobMapKeys.PIPELINE_RUN_OUTPUTS_FILE_SIZE, outputFileSizes);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo
    return StepResult.getStepResultSuccess();
  }
}
