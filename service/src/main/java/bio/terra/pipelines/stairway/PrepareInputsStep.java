package bio.terra.pipelines.stairway;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.Map;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.retry.RetryException;

public class PrepareInputsStep implements Step {
  private final PipelinesService pipelinesService;
  private final Logger logger = LoggerFactory.getLogger(PrepareInputsStep.class);

  public PrepareInputsStep(PipelinesService pipelinesService) {
    this.pipelinesService = pipelinesService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext)
      throws InterruptedException, RetryException {
    // validate and extract parameters from input map
    var inputParameters = flightContext.getInputParameters();
    var workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(
        inputParameters,
        JobMapKeys.PIPELINE_NAME.getKeyName(),
        RunImputationJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS);

    PipelinesEnum pipelineEnum =
        PipelinesEnum.valueOf(
            inputParameters.get(JobMapKeys.PIPELINE_NAME.getKeyName(), String.class));
    Map<String, Object> userProvidedPipelineInputs =
        inputParameters.get(
            RunImputationJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS, new TypeReference<>() {});

    Map<String, Object> allPipelineInputs =
        pipelinesService.constructInputs(pipelineEnum, userProvidedPipelineInputs);

    workingMap.put(RunImputationJobFlightMapKeys.ALL_PIPELINE_INPUTS, allPipelineInputs);
    logger.info("Constructed pipeline inputs: {}", allPipelineInputs);

    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) throws InterruptedException {
    // no undo for this step
    return StepResult.getStepResultSuccess();
  }
}