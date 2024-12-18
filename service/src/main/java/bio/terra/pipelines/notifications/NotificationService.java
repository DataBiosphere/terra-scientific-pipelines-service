package bio.terra.pipelines.notifications;

import static bio.terra.pipelines.app.controller.JobApiUtils.buildApiErrorReport;

import bio.terra.pipelines.app.configuration.internal.NotificationConfiguration;
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
  private final ObjectMapper objectMapper;

  public NotificationService(
      PipelineRunsService pipelineRunsService,
      PipelinesService pipelinesService,
      QuotasService quotasService,
      PubsubService pubsubService,
      NotificationConfiguration notificationConfiguration,
      ObjectMapper objectMapper) {
    this.pipelineRunsService = pipelineRunsService;
    this.pipelinesService = pipelinesService;
    this.quotasService = quotasService;
    this.pubsubService = pubsubService;
    this.notificationConfiguration = notificationConfiguration;
    this.objectMapper = objectMapper;
  }

  /**
   * Pull together the common fields for a notification.
   *
   * @param jobId the job id
   * @param userId the user id
   * @return the base notification object
   */
  public BaseTeaspoonsJobNotification createBaseTeaspoonsJobNotification(
      UUID jobId, String userId) {
    PipelineRun pipelineRun = pipelineRunsService.getPipelineRun(jobId, userId);
    Pipeline pipeline = pipelinesService.getPipelineById(pipelineRun.getPipelineId());
    String pipelineDisplayName = pipeline.getDisplayName();

    // if flight fails before quota steps on user's first run, there won't be a row for them yet
    // in the quotas table
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(userId, pipeline.getName());
    String quotaRemaining = String.valueOf(userQuota.getQuota() - userQuota.getQuotaConsumed());

    return new BaseTeaspoonsJobNotification(
        pipelineRun.getUserId(),
        pipelineDisplayName,
        pipelineRun.getJobId().toString(),
        formatInstantToReadableString(pipelineRun.getCreated()),
        formatInstantToReadableString(pipelineRun.getUpdated()),
        quotaRemaining,
        pipelineRun.getDescription());
  }

  /**
   * Format an Instant as a date time string in UTC using the RFC-1123 date-time formatter, such as
   * 'Tue, 3 Jun 2008 11:05:30 GMT'.
   *
   * @param dateTime the Instant to format
   * @return the formatted date time string
   */
  public String formatInstantToReadableString(Instant dateTime) {
    return dateTime.atZone(ZoneId.of("UTC")).format(DateTimeFormatter.RFC_1123_DATE_TIME);
  }

  /**
   * Create a notification object for a failed job.
   *
   * @param jobId the job id
   * @param userId the user id
   * @param context the flight context
   * @return the notification object
   */
  public TeaspoonsJobFailedNotification createTeaspoonsJobFailedNotification(
      UUID jobId, String userId, FlightContext context) {
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
        createBaseTeaspoonsJobNotification(jobId, userId), errorMessage);
  }

  /**
   * Create a notification object for a successful job.
   *
   * @param jobId the job id
   * @param userId the user id
   * @return the notification object
   */
  public TeaspoonsJobSucceededNotification createTeaspoonsJobSucceededNotification(
      UUID jobId, String userId) {
    String quotaConsumedByJob =
        pipelineRunsService.getPipelineRun(jobId, userId).getQuotaConsumed().toString();
    return new TeaspoonsJobSucceededNotification(
        createBaseTeaspoonsJobNotification(jobId, userId), quotaConsumedByJob);
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
          objectMapper.writeValueAsString(createTeaspoonsJobSucceededNotification(jobId, userId)));
    } catch (IOException | InterruptedException e) {
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
              createTeaspoonsJobFailedNotification(jobId, userId, context)));
    } catch (IOException | InterruptedException e) {
      logger.error("Error sending pipelineRunFailed notification", e);
    }
  }
}
