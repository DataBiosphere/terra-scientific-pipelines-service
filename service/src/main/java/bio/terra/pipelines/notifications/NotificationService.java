package bio.terra.pipelines.notifications;

import static bio.terra.pipelines.app.controller.JobApiUtils.buildApiErrorReport;

import bio.terra.pipelines.app.configuration.internal.NotificationConfiguration;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.generated.model.ApiErrorReport;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.stairway.FlightContext;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.time.Instant;
import java.time.ZoneId;
import java.time.format.DateTimeFormatter;
import java.util.Optional;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Service;

/**
 * Service to encapsulate the logic for composing and sending email notifications to users about
 * completed pipeline runs. Works with the Terra Thurloe service via PubSub messages.
 */
@Service
public class NotificationService {
  private static final Logger logger = LoggerFactory.getLogger(NotificationService.class);

  private final PipelineRunsService pipelineRunsService;
  private final PipelinesService pipelinesService;
  private final QuotasService quotasService;
  private final PubsubService pubsubService;
  private final NotificationConfiguration notificationConfiguration;
  private final PipelineConfigurations pipelineConfigurations;
  private final ObjectMapper objectMapper;

  public NotificationService(
      PipelineRunsService pipelineRunsService,
      PipelinesService pipelinesService,
      QuotasService quotasService,
      PubsubService pubsubService,
      NotificationConfiguration notificationConfiguration,
      PipelineConfigurations pipelineConfigurations,
      ObjectMapper objectMapper) {
    this.pipelineRunsService = pipelineRunsService;
    this.pipelinesService = pipelinesService;
    this.quotasService = quotasService;
    this.pubsubService = pubsubService;
    this.notificationConfiguration = notificationConfiguration;
    this.pipelineConfigurations = pipelineConfigurations;
    this.objectMapper = objectMapper;
  }

  /**
   * Pull together the common fields for a notification.
   *
   * @param jobId the job id
   * @param userId the user id
   * @param context the flight context (only needed for failed notifications)
   * @param isSuccess whether the notification is for a succeeded job; if false, creates a failed
   *     notification
   * @return the base notification object
   */
  private BaseTeaspoonsJobNotification createTeaspoonsJobNotification(
      UUID jobId, String userId, FlightContext context, boolean isSuccess) {
    PipelineRun pipelineRun = pipelineRunsService.getPipelineRun(jobId, userId);
    Pipeline pipeline = pipelinesService.getPipelineById(pipelineRun.getPipelineId());
    String pipelineDisplayName = pipeline.getDisplayName();

    // if flight fails before quota steps on user's first run, there won't be a row for them yet
    // in the quotas table
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(userId, pipeline.getName());
    String quotaRemaining = String.valueOf(userQuota.getQuota() - userQuota.getQuotaConsumed());

    if (isSuccess) { // succeeded
      return new TeaspoonsJobSucceededNotification(
          userId,
          pipelineDisplayName,
          jobId.toString(),
          formatInstantToReadableString(pipelineRun.getCreated()),
          formatInstantToReadableString(pipelineRun.getUpdated()),
          pipelineRun.getQuotaConsumed().toString(),
          quotaRemaining,
          pipelineRun.getDescription(),
          pipelineConfigurations.getCommon().getUserDataTtlDays().toString());
    } else { // failed
      // get exception
      Optional<Exception> exception = context.getResult().getException();
      String errorMessage;
      if (exception.isPresent()) {
        ApiErrorReport errorReport =
            buildApiErrorReport(exception.get()); // use same logic that the status endpoint uses
        errorMessage = errorReport.getMessage();
      } else {
        logger.error(
            "No exception found in flight result for flight {} with status {}",
            context.getFlightId(),
            context.getFlightStatus());
        errorMessage = "Unknown error";
      }
      return new TeaspoonsJobFailedNotification(
          userId,
          pipelineDisplayName,
          jobId.toString(),
          errorMessage,
          formatInstantToReadableString(pipelineRun.getCreated()),
          formatInstantToReadableString(pipelineRun.getUpdated()),
          quotaRemaining,
          pipelineRun.getDescription());
    }
  }

  /**
   * Format an Instant as a date time string in UTC using the RFC-1123 date-time formatter, such as
   * 'Tue, 3 Jun 2008 11:05:30 GMT'.
   *
   * @param dateTime the Instant to format
   * @return the formatted date time string
   */
  protected String formatInstantToReadableString(Instant dateTime) {
    return dateTime.atZone(ZoneId.of("UTC")).format(DateTimeFormatter.RFC_1123_DATE_TIME);
  }

  /**
   * Configure and send a notification that a job has succeeded.
   *
   * @param jobId the job id
   * @param userId the user id
   */
  public void configureAndSendPipelineRunSucceededNotification(UUID jobId, String userId) {
    try {
      pubsubService.publishMessage(
          notificationConfiguration.projectId(),
          notificationConfiguration.topicId(),
          objectMapper.writeValueAsString(
              createTeaspoonsJobNotification(jobId, userId, null, true)));
    } catch (IOException e) {
      logger.error("Error sending pipelineRunSucceeded notification", e);
    }
  }

  /**
   * Configure and send a notification that a job has failed.
   *
   * @param jobId the job id
   * @param userId the user id
   * @param context the flight context
   */
  public void configureAndSendPipelineRunFailedNotification(
      UUID jobId, String userId, FlightContext context) {
    try {
      pubsubService.publishMessage(
          notificationConfiguration.projectId(),
          notificationConfiguration.topicId(),
          objectMapper.writeValueAsString(
              createTeaspoonsJobNotification(jobId, userId, context, false)));
    } catch (IOException e) {
      logger.error("Error sending pipelineRunFailed notification", e);
    }
  }

  /** Configure and send an email notification that a user's quota has changed. */
  public void configureAndSendUserQuotaChangedNotification(
      String userId,
      String pipelineDisplayName,
      int previousQuotaLimit,
      int newQuotaLimit,
      int quotaConsumedByUser,
      int quotaAvailable) {
    TeaspoonsUserQuotaChangedNotification notificationObj =
        new TeaspoonsUserQuotaChangedNotification(
            userId,
            pipelineDisplayName,
            String.valueOf(previousQuotaLimit),
            String.valueOf(newQuotaLimit),
            String.valueOf(quotaConsumedByUser),
            String.valueOf(quotaAvailable));

    try {
      pubsubService.publishMessage(
          notificationConfiguration.projectId(),
          notificationConfiguration.topicId(),
          objectMapper.writeValueAsString(notificationObj));
    } catch (IOException e) {
      logger.error("Error sending user quota changed notification to pubsub", e);
    }
  }
}
