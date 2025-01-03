package bio.terra.pipelines.notifications;

import lombok.Getter;

/** Base class for Teaspoons job notifications. Contains common fields for all job notifications. */
@Getter
public abstract class BaseTeaspoonsJobNotification {
  protected String notificationType;
  protected String recipientUserId;
  protected String pipelineDisplayName;
  protected String jobId;
  protected String timeSubmitted;
  protected String timeCompleted;
  protected String quotaRemaining;
  protected String quotaConsumedByJob;
  protected String userDescription;

  protected BaseTeaspoonsJobNotification(
      String recipientUserId,
      String pipelineDisplayName,
      String jobId,
      String timeSubmitted,
      String timeCompleted,
      String quotaRemaining,
      String userDescription) {
    this.recipientUserId = recipientUserId;
    this.pipelineDisplayName = pipelineDisplayName;
    this.jobId = jobId;
    this.timeSubmitted = timeSubmitted;
    this.timeCompleted = timeCompleted;
    this.quotaRemaining = quotaRemaining;
    this.userDescription = userDescription;
  }
}
