package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.common.utils.FlightUtils.flightMapKeyIsTrue;

import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.notifications.NotificationService;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightStatus;
import bio.terra.stairway.HookAction;
import bio.terra.stairway.StairwayHook;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

/**
 * A {@link StairwayHook} that sends a Job Failed Notification email via pubsub/Thurloe upon flight
 * failure.
 *
 * <p>This hook action will only run if the flight's input parameters contain the JobMapKeys key for
 * DO_SEND_JOB_FAILURE_NOTIFICATION_HOOK and the flight's status is not SUCCESS.
 *
 * <p>The JobMapKeys key for PIPELINE_NAME is required to send the notification.
 */
@Component
public class StairwaySendFailedJobNotificationHook implements StairwayHook {
  private final NotificationService notificationService;
  private static final Logger logger =
      LoggerFactory.getLogger(StairwaySendFailedJobNotificationHook.class);

  public StairwaySendFailedJobNotificationHook(NotificationService notificationService) {
    this.notificationService = notificationService;
  }

  @Override
  public HookAction endFlight(FlightContext context) {

    FlightMap inputParameters = context.getInputParameters();

    if (flightMapKeyIsTrue(inputParameters, JobMapKeys.DO_SEND_JOB_FAILURE_NOTIFICATION_HOOK)
        && context.getFlightStatus() != FlightStatus.SUCCESS) {
      logger.info(
          "Flight has status {}, sending failed job notification email", context.getFlightStatus());

      FlightUtils.validateRequiredEntries(inputParameters, JobMapKeys.USER_ID);

      UUID jobId = UUID.fromString(context.getFlightId());
      String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);

      // send email notification
      notificationService.configureAndSendPipelineRunFailedNotification(jobId, userId, context);
    }
    return HookAction.CONTINUE;
  }
}
