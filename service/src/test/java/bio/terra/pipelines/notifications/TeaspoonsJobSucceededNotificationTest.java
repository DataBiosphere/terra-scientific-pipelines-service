package bio.terra.pipelines.notifications;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.testutils.BaseTest;
import org.junit.jupiter.api.Test;

class TeaspoonsJobSucceededNotificationTest extends BaseTest {

  String recipientUserId = "test-recipient-user-id";
  String pipelineDisplayName = "test-pipeline-display-name";
  String jobId = "test-job-id";
  String timeSubmitted = "test-time-submitted";
  String timeCompleted = "test-time-completed";
  String quotaConsumedByJob = "test-quota-consumed-by-job";
  String quotaRemaining = "test-quota-remaining";
  String userDescription = "test-user-description";
  String userDataTtl = "test-user-data-ttl";

  @Test
  void teaspoonsJobSucceededNotificationFullConstructor() {
    TeaspoonsJobSucceededNotification notification =
        new TeaspoonsJobSucceededNotification(
            recipientUserId,
            pipelineDisplayName,
            jobId,
            timeSubmitted,
            timeCompleted,
            quotaConsumedByJob,
            quotaRemaining,
            userDescription,
            userDataTtl);
    assertEquals("TeaspoonsJobSucceededNotification", notification.getNotificationType());
    assertEquals(recipientUserId, notification.getRecipientUserId());
    assertEquals(pipelineDisplayName, notification.getPipelineDisplayName());
    assertEquals(jobId, notification.getJobId());
    assertEquals(timeSubmitted, notification.getTimeSubmitted());
    assertEquals(timeCompleted, notification.getTimeCompleted());
    assertEquals(quotaConsumedByJob, notification.getQuotaConsumedByJob());
    assertEquals(quotaRemaining, notification.getQuotaRemaining());
    assertEquals(userDescription, notification.getUserDescription());
    assertEquals(userDataTtl, notification.getUserDataTtl());
  }
}
