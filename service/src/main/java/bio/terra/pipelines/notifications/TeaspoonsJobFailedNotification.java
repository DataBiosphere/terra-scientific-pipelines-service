package bio.terra.pipelines.notifications;

import lombok.Getter;

@Getter
public class TeaspoonsJobFailedNotification extends BaseTeaspoonsJobNotification {
  private static final String NOTIFICATION_TYPE = "TeaspoonsJobFailedNotification";
  private static final String QUOTA_CONSUMED_BY_FAILED_JOB = "0";
  public final String notificationType;
  public final String errorMessage;
  public final String quotaConsumedByJob;

  public TeaspoonsJobFailedNotification(
      BaseTeaspoonsJobNotification baseTeaspoonsJobNotification, String errorMessage) {
    super(
        baseTeaspoonsJobNotification.recipientUserId,
        baseTeaspoonsJobNotification.pipelineDisplayName,
        baseTeaspoonsJobNotification.jobId,
        baseTeaspoonsJobNotification.timeSubmitted,
        baseTeaspoonsJobNotification.timeCompleted,
        baseTeaspoonsJobNotification.quotaRemaining,
        baseTeaspoonsJobNotification.userDescription);
    this.notificationType = NOTIFICATION_TYPE;
    this.errorMessage = errorMessage;
    this.quotaConsumedByJob = QUOTA_CONSUMED_BY_FAILED_JOB;
  }

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
    this.errorMessage = errorMessage;
    this.quotaConsumedByJob = QUOTA_CONSUMED_BY_FAILED_JOB;
  }
}
