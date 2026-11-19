package bio.terra.pipelines.notifications.model;

import lombok.Getter;

@Getter
public class TeaspoonsUserQuotaChangedNotification extends BaseTeaspoonsNotification {
  private final String recipientUserId;
  private final String customMessage;
  private final String previousQuotaAllocated;
  private final String updatedQuotaAllocated;
  private final String previousQuotaConsumed;
  private final String updatedQuotaConsumed;
  private final String quotaAvailable;

  public TeaspoonsUserQuotaChangedNotification(
      String recipientUserId,
      String pipelineDisplayName,
      String customMessage,
      String previousQuotaAllocated,
      String updatedQuotaAllocated,
      String previousQuotaConsumed,
      String updatedQuotaConsumed,
      String quotaAvailable) {
    this.notificationType = "TeaspoonsUserQuotaChangedNotification";
    this.recipientUserId = recipientUserId;
    this.pipelineDisplayName = pipelineDisplayName;
    this.customMessage = customMessage;
    this.previousQuotaAllocated = previousQuotaAllocated;
    this.updatedQuotaAllocated = updatedQuotaAllocated;
    this.previousQuotaConsumed = previousQuotaConsumed;
    this.updatedQuotaConsumed = updatedQuotaConsumed;
    this.quotaAvailable = quotaAvailable;
  }
}
