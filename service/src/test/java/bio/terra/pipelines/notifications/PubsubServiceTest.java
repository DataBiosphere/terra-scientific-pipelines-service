package bio.terra.pipelines.notifications;

import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.io.IOException;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class PubsubServiceTest extends BaseEmbeddedDbTest {
  @Autowired PubsubService pubsubService;

  @Test
  void publishMessage() throws IOException, InterruptedException {
    // setup
    String message = "test message";
    String topicId = "testTopicId";
    String projectId = "testProjectId";

    // lol this passes the test with no mocking since we catch the exception
    pubsubService.publishMessage(projectId, topicId, message);
  }
}
