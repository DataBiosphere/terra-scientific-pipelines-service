package bio.terra.pipelines.notifications;

public record TeaspoonsJobSucceededNotification(
    String notificationType,
    String recipientUserId,
    String pipelineDisplayName,
    String jobId,
    String timeSubmitted,
    String timeCompleted,
    String quotaConsumedByJob,
    String quotaRemaining,
    String userDescription) {

  public TeaspoonsJobSucceededNotification(
      String recipientUserId,
      String pipelineDisplayName,
      String jobId,
      String timeSubmitted,
      String timeCompleted,
      String quotaConsumedByJob,
      String quotaRemaining,
      String userDescription) {
    this(
        "TeaspoonsJobSucceededNotification",
        recipientUserId,
        pipelineDisplayName,
        jobId,
        timeSubmitted,
        timeCompleted,
        quotaConsumedByJob,
        quotaRemaining,
        userDescription);
  }
}
