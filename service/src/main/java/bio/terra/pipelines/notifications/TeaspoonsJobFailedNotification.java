package bio.terra.pipelines.notifications;

public record TeaspoonsJobFailedNotification(
    String notificationType,
    String recipientUserId,
    String pipelineDisplayName,
    String jobId,
    String errorMessage,
    String timeSubmitted,
    String timeCompleted,
    String quotaConsumedByJob,
    String quotaRemaining,
    String userDescription) {

  public TeaspoonsJobFailedNotification(
      String recipientUserId,
      String pipelineDisplayName,
      String jobId,
      String errorMessage,
      String timeSubmitted,
      String timeCompleted,
      String quotaConsumedByJob,
      String quotaRemaining,
      String userDescription) {
    this(
        "TeaspoonsJobFailedNotification",
        recipientUserId,
        pipelineDisplayName,
        jobId,
        errorMessage,
        timeSubmitted,
        timeCompleted,
        quotaConsumedByJob,
        quotaRemaining,
        userDescription);
  }
}
