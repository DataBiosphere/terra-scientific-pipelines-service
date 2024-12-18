package bio.terra.pipelines.notifications;

import lombok.Getter;

@Getter
public class TeaspoonsJobSucceededNotification extends BaseTeaspoonsJobNotification {
  private static final String NOTIFICATION_TYPE = "TeaspoonsJobSucceededNotification";
  public final String notificationType;
  public final String quotaConsumedByJob;

  public TeaspoonsJobSucceededNotification(
      BaseTeaspoonsJobNotification baseTeaspoonsJobNotification, String quotaConsumedByJob) {
    super(
        baseTeaspoonsJobNotification.recipientUserId,
        baseTeaspoonsJobNotification.pipelineDisplayName,
        baseTeaspoonsJobNotification.jobId,
        baseTeaspoonsJobNotification.timeSubmitted,
        baseTeaspoonsJobNotification.timeCompleted,
        baseTeaspoonsJobNotification.quotaRemaining,
        baseTeaspoonsJobNotification.userDescription);
    this.notificationType = NOTIFICATION_TYPE;
    this.quotaConsumedByJob = quotaConsumedByJob;
  }

  public TeaspoonsJobSucceededNotification(
      String recipientUserId,
      String pipelineDisplayName,
      String jobId,
      String timeSubmitted,
      String timeCompleted,
      String quotaConsumedByJob,
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
    this.quotaConsumedByJob = quotaConsumedByJob;
  }
}
