package bio.terra.pipelines.notifications.model;

import lombok.Getter;

@Getter
@SuppressWarnings({"java:S107"}) // Disable "Methods should not have too many parameters"
public class TeaspoonsJobFailedNotification extends TeaspoonsJobNotification {
  private static final String NOTIFICATION_TYPE = "TeaspoonsJobFailedNotification";
  private static final String QUOTA_CONSUMED_BY_FAILED_JOB = "0";
  private final String errorMessage;

  public TeaspoonsJobFailedNotification(
      String recipientUserId,
      String pipelineDisplayName,
      String jobId,
      String errorMessage,
      String timeSubmitted,
      String timeCompleted,
      String quotaRemaining,
      String userDescription) {
    super(
        recipientUserId,
        pipelineDisplayName,
        jobId,
        timeSubmitted,
        timeCompleted,
        quotaRemaining,
        userDescription);
    this.notificationType = NOTIFICATION_TYPE;
    this.quotaConsumedByJob = QUOTA_CONSUMED_BY_FAILED_JOB;
    this.errorMessage = errorMessage;
  }
}
