package bio.terra.pipelines.notifications;

import com.google.api.core.ApiFuture;
import com.google.api.core.ApiFutureCallback;
import com.google.api.core.ApiFutures;
import com.google.api.gax.rpc.ApiException;
import com.google.cloud.pubsub.v1.Publisher;
import com.google.common.util.concurrent.MoreExecutors;
import com.google.protobuf.ByteString;
import com.google.pubsub.v1.PubsubMessage;
import com.google.pubsub.v1.TopicName;
import java.io.IOException;
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

    // Add an asynchronous callback to handle success / failure
    // Registers a callback to be run when the ApiFuture's computation is complete or, if the
    // computation is already complete, immediately.
    ApiFutures.addCallback(
        future, handleCompletedPubsubMessagePublishing(logger), MoreExecutors.directExecutor());
  }

  ApiFutureCallback<String> handleCompletedPubsubMessagePublishing(Logger logger) {
    return new ApiFutureCallback<String>() {

      @Override
      public void onFailure(Throwable throwable) {
        String errorMessage = "";
        if (throwable instanceof ApiException apiException) {
          // details on the API exception
          errorMessage =
              "Google API exception: status code %s, is retryable: %s"
                  .formatted(apiException.getStatusCode().getCode(), apiException.isRetryable());
        }
        logger.error(
            "Error publishing message to Google PubSub: {}; {}",
            errorMessage,
            throwable.getMessage());
      }

      @Override
      public void onSuccess(String messageId) {
        // Once published, returns server-assigned message ids (unique within the topic)
        logger.info("Published message ID: {}", messageId);
      }
    };
  }
}
