package bio.terra.pipelines.notifications;

import com.google.api.core.ApiFuture;
import com.google.api.gax.rpc.ApiException;
import com.google.cloud.pubsub.v1.Publisher;
import com.google.protobuf.ByteString;
import com.google.pubsub.v1.PubsubMessage;
import com.google.pubsub.v1.TopicName;
import java.io.IOException;
import java.util.concurrent.TimeUnit;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Service;

/** Service to interact with Pubsub. Used by NotificationService. */
@Service
public class PubsubService {
  private static final Logger logger = LoggerFactory.getLogger(PubsubService.class);

  protected static Publisher publisher;

  /**
   * Initialize a publisher for a Google PubSub topic. Does nothing if the publisher already exists.
   */
  protected static void initPublisher(TopicName topicName) throws IOException {
    if (publisher != null) {
      logger.info("Publisher already exists for Google PubSub topicName {}", topicName);
      return;
    }
    try {
      logger.info("Creating publisher for Google PubSub topicName {}", topicName);
      // Create a publisher instance with default settings bound to the topic
      publisher = Publisher.newBuilder(topicName).build();
    } catch (IOException e) {
      logger.error("Error creating publisher for Google PubSub topicName {}", topicName);
      throw e;
    }
  }

  /** Publish a message to a Google PubSub topic. */
  public void publishMessage(String projectId, String topicId, String message) throws IOException {
    TopicName topicName = TopicName.of(projectId, topicId);
    logger.info(
        "Publishing message to Google PubSub projectId {}, topicId {}: {}",
        projectId,
        topicId,
        message);

    initPublisher(topicName);

    ByteString data = ByteString.copyFromUtf8(message);
    PubsubMessage pubsubMessage = PubsubMessage.newBuilder().setData(data).build();

    // Once published, returns a server-assigned message id (unique within the topic)
    ApiFuture<String> future = publisher.publish(pubsubMessage);

    try {
      // Wait on any pending publish requests with a 30-second timeout
      String messageId = future.get(30, TimeUnit.SECONDS);
      logger.info("Published message ID: {}", messageId);
    } catch (Exception e) {
      String errorMessage;
      if (e instanceof ApiException apiException) {
        // details on the API exception
        errorMessage =
            "Google API exception: status code %s, is retryable: %s"
                .formatted(apiException.getStatusCode().getCode(), apiException.isRetryable());
        logger.error(
            "Error publishing message to Google PubSub: {}; {}", errorMessage, e.getMessage());
      } else if (e instanceof InterruptedException) {
        errorMessage = "Thread was interrupted";
        logger.error(
            "Error publishing message to Google PubSub: {}; {}", errorMessage, e.getMessage());
        Thread.currentThread().interrupt();
      } else {
        logger.error("Error publishing message to Google PubSub: {}", e.getMessage());
      }
    }
  }
}
