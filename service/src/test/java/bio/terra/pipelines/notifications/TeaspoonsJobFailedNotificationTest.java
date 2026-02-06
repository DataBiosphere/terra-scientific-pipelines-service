package bio.terra.pipelines.notifications;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.notifications.model.TeaspoonsJobFailedNotification;
import bio.terra.pipelines.testutils.BaseTest;
import org.junit.jupiter.api.Test;

class TeaspoonsJobFailedNotificationTest extends BaseTest {
  String recipientUserId = "test-recipient-user-id";
  String pipelineDisplayName = "test-pipeline-display-name";
  String jobId = "test-job-id";
  String errorMessage = "test-error-message";
  String timeSubmitted = "test-time-submitted";
  String timeCompleted = "test-time-completed";
  String quotaRemaining = "test-quota-remaining";
  String userDescription = "test-user-description";

  @Test
  void teaspoonsJobFailedNotificationFullConstructor() {

    TeaspoonsJobFailedNotification notification =
        new TeaspoonsJobFailedNotification(
            recipientUserId,
            pipelineDisplayName,
            jobId,
            errorMessage,
            timeSubmitted,
            timeCompleted,
            quotaRemaining,
            userDescription);
    assertEquals("TeaspoonsJobFailedNotification", notification.getNotificationType());
    assertEquals(recipientUserId, notification.getRecipientUserId());
    assertEquals(pipelineDisplayName, notification.getPipelineDisplayName());
    assertEquals(jobId, notification.getJobId());
    assertEquals(errorMessage, notification.getErrorMessage());
    assertEquals(timeSubmitted, notification.getTimeSubmitted());
    assertEquals(timeCompleted, notification.getTimeCompleted());
    assertEquals("0", notification.getQuotaConsumedByJob());
    assertEquals(quotaRemaining, notification.getQuotaRemaining());
    assertEquals(userDescription, notification.getUserDescription());
  }
}
