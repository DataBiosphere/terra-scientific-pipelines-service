package bio.terra.pipelines.notifications.model;

import bio.terra.pipelines.model.UserQuotaTotals;
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
      UserQuotaTotals previousQuotaTotals,
      UserQuotaTotals updatedQuotaTotals,
      String quotaAvailable) {
    this.notificationType = "TeaspoonsUserQuotaChangedNotification";
    this.recipientUserId = recipientUserId;
    this.pipelineDisplayName = pipelineDisplayName;
    this.customMessage = customMessage;
    this.previousQuotaAllocated = String.valueOf(previousQuotaTotals.totalAllocated());
    this.updatedQuotaAllocated = String.valueOf(updatedQuotaTotals.totalAllocated());
    this.previousQuotaConsumed = String.valueOf(previousQuotaTotals.totalConsumed());
    this.updatedQuotaConsumed = String.valueOf(updatedQuotaTotals.totalConsumed());
    this.quotaAvailable = quotaAvailable;
  }
}
