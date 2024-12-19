package bio.terra.pipelines.notifications;

import bio.terra.pipelines.app.configuration.internal.NotificationConfiguration;
import bio.terra.pipelines.db.entities.PipelineRun;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.google.common.annotations.VisibleForTesting;
import jakarta.annotation.PostConstruct;
import java.io.IOException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Service;

/**
 * Service to encapsulate the logic for sending email notifications to users. Works with the Terra
 * Thurloe service via PubSub messages.
 */
@Service
public class NotificationService {
  private static final Logger logger = LoggerFactory.getLogger(NotificationService.class);

  private final PubsubService pubsubService;
  private final NotificationConfiguration notificationConfiguration;
  private final ObjectMapper objectMapper;

  public NotificationService(
      PubsubService pubsubService,
      NotificationConfiguration notificationConfiguration,
      ObjectMapper objectMapper) {
    this.pubsubService = pubsubService;
    this.notificationConfiguration = notificationConfiguration;
    this.objectMapper = objectMapper;
  }

  @VisibleForTesting
  @PostConstruct
  protected void createTopic() {
    try {
      pubsubService.createTopic(
          notificationConfiguration.projectId(), notificationConfiguration.topicId());
    } catch (IOException e) {
      logger.warn("Error creating notification topic", e);
    }
  }

  public void sendPipelineRunSucceededNotification(
      PipelineRun pipelineRun, String pipelineDisplayName, String quotaRemaining) {
    try {
      pubsubService.publishMessage(
          notificationConfiguration.projectId(),
          notificationConfiguration.topicId(),
          objectMapper.writeValueAsString(
              new TeaspoonsJobSucceededNotification(
                  pipelineRun.getUserId(),
                  pipelineDisplayName,
                  pipelineRun.getJobId().toString(),
                  pipelineRun.getCreated().toString(),
                  pipelineRun.getUpdated().toString(),
                  pipelineRun.getQuotaConsumed().toString(),
                  quotaRemaining,
                  pipelineRun.getDescription())));
    } catch (IOException e) {
      logger.error("Error sending pipelineRunSucceeded notification", e);
    }
  }

  public void sendPipelineRunFailedNotification(
      PipelineRun pipelineRun,
      String pipelineDisplayName,
      String quotaRemaining,
      String errorMessage) {
    try {
      pubsubService.publishMessage(
          notificationConfiguration.projectId(),
          notificationConfiguration.topicId(),
          objectMapper.writeValueAsString(
              new TeaspoonsJobFailedNotification(
                  pipelineRun.getUserId(),
                  pipelineDisplayName,
                  pipelineRun.getJobId().toString(),
                  errorMessage,
                  pipelineRun.getCreated().toString(),
                  pipelineRun.getUpdated().toString(),
                  "0",
                  quotaRemaining,
                  pipelineRun.getDescription())));
    } catch (IOException e) {
      logger.error("Error sending pipelineRunFailed notification", e);
    }
  }
}
