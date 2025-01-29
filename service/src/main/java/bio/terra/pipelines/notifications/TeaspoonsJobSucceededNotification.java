package bio.terra.pipelines.notifications;

import lombok.Getter;

@Getter
@SuppressWarnings({"java:S107"}) // Disable "Methods should not have too many parameters"
public class TeaspoonsJobSucceededNotification extends BaseTeaspoonsJobNotification {
  private static final String NOTIFICATION_TYPE = "TeaspoonsJobSucceededNotification";
  private final String userDataTtlDays;

  public TeaspoonsJobSucceededNotification(
      String recipientUserId,
      String pipelineDisplayName,
      String jobId,
      String timeSubmitted,
      String timeCompleted,
      String quotaConsumedByJob,
      String quotaRemaining,
      String userDescription,
      String userDataTtlDays) {
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
    this.userDataTtlDays = userDataTtlDays;
  }
}
