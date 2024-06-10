package bio.terra.pipelines.stairway.imputation;

import bio.terra.cbas.model.*;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.cbas.CbasServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.stairway.*;
import java.util.List;
import java.util.UUID;
import java.util.concurrent.TimeUnit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This step polls cbas for a run set submission until all runs are in a finalized state. If not all
 * runs are in a finalized state, the step will sleep for
 * bio.terra.pipelines.app.configuration.internal.ImputationConfiguration#cromwellSubmissionPollingIntervalInSeconds
 * seconds before polling again. Once all runs are finalized then it will see if they are all
 * successful and if so will succeed otherwise will fail.
 *
 * <p>this step expects cbas uri and run set id to provided in the working map
 */
public class PollCromwellRunSetStatusStep implements Step {
  private final CbasService cbasService;
  private final SamService samService;
  private final ImputationConfiguration imputationConfiguration;
  private final Logger logger = LoggerFactory.getLogger(PollCromwellRunSetStatusStep.class);

  public PollCromwellRunSetStatusStep(
      CbasService cbasService,
      SamService samService,
      ImputationConfiguration imputationConfiguration) {
    this.cbasService = cbasService;
    this.samService = samService;
    this.imputationConfiguration = imputationConfiguration;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) throws InterruptedException {
    // validate and extract parameters from working map
    FlightMap workingMap = flightContext.getWorkingMap();
    FlightUtils.validateRequiredEntries(
        workingMap,
        RunImputationJobFlightMapKeys.CBAS_URI,
        RunImputationJobFlightMapKeys.RUN_SET_ID);

    String cbasUri = workingMap.get(RunImputationJobFlightMapKeys.CBAS_URI, String.class);
    UUID runSetId = workingMap.get(RunImputationJobFlightMapKeys.RUN_SET_ID, UUID.class);

    // poll until all runs are in a finalized state
    RunLogResponse runLogResponse = null;
    boolean stillRunning = true;
    try {
      while (stillRunning) {
        runLogResponse =
            cbasService.getRunsForRunSet(
                cbasUri, samService.getTspsServiceAccountToken(), runSetId);
        stillRunning = CbasService.containsRunningRunLog(runLogResponse);
        if (stillRunning) {
          logger.info(
              "Polling Started, sleeping for {} seconds",
              imputationConfiguration.getCromwellSubmissionPollingIntervalInSeconds());
          TimeUnit.SECONDS.sleep(
              imputationConfiguration.getCromwellSubmissionPollingIntervalInSeconds());
        }
      }
    } catch (CbasServiceApiException e) {
      return new StepResult(StepStatus.STEP_RESULT_FAILURE_RETRY, e);
    }

    // if there are any non-successful logs, fatally fail the step
    List<RunLog> failedRunLogs =
        runLogResponse.getRuns().stream()
            .filter(runLog -> !runLog.getState().equals(RunState.COMPLETE))
            .toList();
    if (failedRunLogs.isEmpty()) {
      return StepResult.getStepResultSuccess();
    } else {
      return new StepResult(
          StepStatus.STEP_RESULT_FAILURE_FATAL,
          new InternalServerErrorException("Not all runs succeeded for run set " + runSetId));
    }
  }

  @Override
  public StepResult undoStep(FlightContext context) {
    // nothing to undo; there's nothing to undo about polling a cromwell run set
    return StepResult.getStepResultSuccess();
  }
}
