package bio.terra.pipelines.notifications;

import lombok.Getter;

@Getter
public class TeaspoonsUserQuotaChangedNotification {

  protected String notificationType;
  protected String recipientUserId;
  protected String pipelineDisplayName;
  protected String previousQuotaLimit;
  protected String newQuotaLimit;
  protected String quotaConsumedByUser;
  protected String quotaAvailable;

  public TeaspoonsUserQuotaChangedNotification(
      String recipientUserId,
      String pipelineDisplayName,
      String previousQuotaLimit,
      String newQuotaLimit,
      String quotaConsumedByUser,
      String quotaAvailable) {
    this.recipientUserId = recipientUserId;
    this.pipelineDisplayName = pipelineDisplayName;
    this.previousQuotaLimit = previousQuotaLimit;
    this.newQuotaLimit = newQuotaLimit;
    this.quotaConsumedByUser = quotaConsumedByUser;
    this.quotaAvailable = quotaAvailable;
  }
}
