package bio.terra.pipelines.notifications.model;

import lombok.Getter;

/**
 * Parent class for Teaspoons job notifications. Contains common fields for all job notifications.
 */
@Getter
public abstract class TeaspoonsJobNotification extends BaseTeaspoonsNotification {
  protected String recipientUserId;
  protected String jobId;
  protected String timeSubmitted;
  protected String timeCompleted;
  protected String quotaRemaining;
  protected String quotaConsumedByJob;
  protected String userDescription;

  protected TeaspoonsJobNotification(
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
