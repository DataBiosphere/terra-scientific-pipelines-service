package bio.terra.pipelines.notifications;

import com.google.api.core.ApiFuture;
import com.google.api.core.ApiFutureCallback;
import com.google.api.core.ApiFutures;
import com.google.api.gax.rpc.ApiException;
import com.google.cloud.pubsub.v1.Publisher;
import com.google.cloud.pubsub.v1.TopicAdminClient;
import com.google.common.util.concurrent.MoreExecutors;
import com.google.protobuf.ByteString;
import com.google.pubsub.v1.PubsubMessage;
import com.google.pubsub.v1.Topic;
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

  public void createTopic(String projectId, String topicId) throws IOException {
    try (TopicAdminClient topicAdminClient = TopicAdminClient.create()) {
      TopicName topicName = TopicName.of(projectId, topicId);
      if (topicAdminClient.getTopic(topicName) == null) {
        Topic topic = topicAdminClient.createTopic(topicName);
        logger.info("Created topic: {}", topic.getName());
      }
    } catch (IOException e) {
      logger.error("Error creating topic", e);
    }
  }

  public void publishMessage(String projectId, String topicId, String message)
      throws IOException, InterruptedException {
    //    TopicName topicName = TopicName.of(projectId, topicId);
    //    var publisher = Publisher.newBuilder(topicName).build();
    //
    // publisher.publish(PubsubMessage.newBuilder().setData(ByteString.copyFromUtf8(message)).build());
    //    publisher.awaitTermination()

    TopicName topicName = TopicName.of(projectId, topicId);
    Publisher publisher = null;

    try {
      // Create a publisher instance with default settings bound to the topic
      publisher = Publisher.newBuilder(topicName).build();

      ByteString data = ByteString.copyFromUtf8(message);
      PubsubMessage pubsubMessage = PubsubMessage.newBuilder().setData(data).build();

      // Once published, returns a server-assigned message id (unique within the topic)
      ApiFuture<String> future = publisher.publish(pubsubMessage);

      // Add an asynchronous callback to handle success / failure
      ApiFutures.addCallback(
          future,
          new ApiFutureCallback<String>() {

            @Override
            public void onFailure(Throwable throwable) {
              if (throwable instanceof ApiException) {
                ApiException apiException = ((ApiException) throwable);
                // details on the API exception
                logger.warn(
                    "Google API exception! status code: {}, is retryable: {}",
                    apiException.getStatusCode().getCode(),
                    apiException.isRetryable());
              }
              logger.error("Error publishing message to Google PubSub: {}", message);
            }

            @Override
            public void onSuccess(String messageId) {
              // Once published, returns server-assigned message ids (unique within the topic)
              logger.info("Published message ID: {}", messageId);
            }
          },
          MoreExecutors.directExecutor());

    } finally {
      if (publisher != null) {
        // When finished with the publisher, shutdown to free up resources.
        publisher.shutdown();
        publisher.awaitTermination(1, TimeUnit.MINUTES);
      }
    }
  }
}
