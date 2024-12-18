package bio.terra.pipelines.notifications;

import lombok.Getter;

@Getter
public class BaseTeaspoonsJobNotification {
  public final String recipientUserId;
  public final String pipelineDisplayName;
  public final String jobId;
  public final String timeSubmitted;
  public final String timeCompleted;
  public final String quotaRemaining;
  public final String userDescription;

  public BaseTeaspoonsJobNotification(
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
