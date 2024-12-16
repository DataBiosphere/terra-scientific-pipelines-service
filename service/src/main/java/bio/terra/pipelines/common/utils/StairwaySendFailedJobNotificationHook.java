package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.app.controller.JobApiUtils.buildApiErrorReport;
import static bio.terra.pipelines.common.utils.FlightUtils.flightMapKeyIsTrue;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.generated.model.ApiErrorReport;
import bio.terra.pipelines.notifications.NotificationService;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightStatus;
import bio.terra.stairway.HookAction;
import bio.terra.stairway.StairwayHook;
import java.util.Optional;
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
  private final PipelineRunsService pipelineRunsService;
  private final PipelinesService pipelinesService;
  private final QuotasService quotasService;
  private final NotificationService notificationService;
  private static final Logger logger =
      LoggerFactory.getLogger(StairwaySendFailedJobNotificationHook.class);

  public StairwaySendFailedJobNotificationHook(
      PipelineRunsService pipelineRunsService,
      NotificationService notificationService,
      PipelinesService pipelinesService,
      QuotasService quotasService) {
    this.pipelineRunsService = pipelineRunsService;
    this.notificationService = notificationService;
    this.pipelinesService = pipelinesService;
    this.quotasService = quotasService;
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

      PipelineRun pipelineRun = pipelineRunsService.getPipelineRun(jobId, userId);
      Pipeline pipeline = pipelinesService.getPipelineById(pipelineRun.getPipelineId());
      String pipelineDisplayName = pipeline.getDisplayName();
      // if flight fails before quota steps on user's first run, there won't be a row for them yet
      // in the quotas table
      UserQuota userQuota =
          quotasService.getOrCreateQuotaForUserAndPipeline(userId, pipeline.getName());
      String quotaRemaining = String.valueOf(userQuota.getQuota() - userQuota.getQuotaConsumed());

      // get exception
      Optional<Exception> exception = context.getResult().getException();
      String errorMessage;
      if (exception.isPresent()) {
        ApiErrorReport errorReport =
            buildApiErrorReport(exception.get()); // use same logic that the status endpoint uses
        errorMessage = errorReport.getMessage();
      } else {
        // should this just fail?
        logger.warn(
            "No exception found in flight result for flight {} with status {}",
            context.getFlightId(),
            context.getFlightStatus());
        errorMessage = "Unknown error";
      }

      // send email notification
      notificationService.sendPipelineRunFailedNotification(
          pipelineRun, pipelineDisplayName, quotaRemaining, errorMessage);
    }
    return HookAction.CONTINUE;
  }
}
