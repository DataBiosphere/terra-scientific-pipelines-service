package bio.terra.pipelines.stairway.steps.common;

import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.notifications.NotificationService;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.Step;
import bio.terra.stairway.StepResult;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Step to send an email notification that a job has succeeded. This step cannot fail (we catch all
 * exceptions).
 */
public class SendJobSucceededNotificationStep implements Step {
  private final NotificationService notificationService;
  private final Logger logger = LoggerFactory.getLogger(CompletePipelineRunStep.class);

  public SendJobSucceededNotificationStep(NotificationService notificationService) {
    this.notificationService = notificationService;
  }

  @Override
  public StepResult doStep(FlightContext flightContext) {
    // we place the entire logic of this step in a try-catch so that it cannot fail
    try {
      // validate and extract parameters from input map
      var inputParameters = flightContext.getInputParameters();
      FlightUtils.validateRequiredEntries(inputParameters, JobMapKeys.USER_ID);

      UUID jobId = UUID.fromString(flightContext.getFlightId());
      String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);

      // send email notification
      notificationService.configureAndSendPipelineRunSucceededNotification(jobId, userId);
    } catch (Exception e) {
      logger.error("Failed to send email notification", e);
    }
    return StepResult.getStepResultSuccess();
  }

  @Override
  public StepResult undoStep(FlightContext flightContext) {
    // nothing to undo
    return StepResult.getStepResultSuccess();
  }
}
