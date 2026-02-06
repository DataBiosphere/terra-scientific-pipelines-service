package bio.terra.pipelines.notifications.model;

import lombok.Getter;

@Getter
public class TeaspoonsUserQuotaChangedNotification extends BaseTeaspoonsNotification {
  private final String recipientUserId;
  private final String previousQuotaLimit;
  private final String newQuotaLimit;
  private final String quotaConsumedByUser;
  private final String quotaAvailable;

  public TeaspoonsUserQuotaChangedNotification(
      String recipientUserId,
      String pipelineDisplayName,
      String previousQuotaLimit,
      String newQuotaLimit,
      String quotaConsumedByUser,
      String quotaAvailable) {
    this.notificationType = "TeaspoonsUserQuotaChangedNotification";
    this.recipientUserId = recipientUserId;
    this.pipelineDisplayName = pipelineDisplayName;
    this.previousQuotaLimit = previousQuotaLimit;
    this.newQuotaLimit = newQuotaLimit;
    this.quotaConsumedByUser = quotaConsumedByUser;
    this.quotaAvailable = quotaAvailable;
  }
}
