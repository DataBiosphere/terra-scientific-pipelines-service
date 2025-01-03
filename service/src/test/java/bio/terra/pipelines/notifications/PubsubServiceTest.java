package bio.terra.pipelines.notifications;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.mockStatic;
import static org.mockito.Mockito.times;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import com.google.cloud.pubsub.v1.Publisher;
import com.google.pubsub.v1.TopicName;
import java.io.IOException;
import org.junit.jupiter.api.Test;
import org.mockito.MockedStatic;
import org.springframework.beans.factory.annotation.Autowired;

class PubsubServiceTest extends BaseEmbeddedDbTest {
  @Autowired PubsubService pubsubService;

  @Test
  void initPublisher() throws IOException {
    TopicName topicName = TopicName.of("projectId", "topicId");
    try (MockedStatic<Publisher> mockedStaticPublisher = mockStatic(Publisher.class)) {
      Publisher.Builder mockBuilder = mock(Publisher.Builder.class);
      mockedStaticPublisher.when(() -> Publisher.newBuilder(topicName)).thenReturn(mockBuilder);
      Publisher mockPublisher = mock(Publisher.class);
      when(mockBuilder.build()).thenReturn(mockPublisher);

      // to start, publisher should be null
      assertNull(PubsubService.publisher);

      // calling init once should set the publisher
      PubsubService.initPublisher(topicName);
      assertEquals(mockPublisher, PubsubService.publisher);

      // calling init again should not call build again
      PubsubService.initPublisher(topicName);
      assertEquals(mockPublisher, PubsubService.publisher);
      verify(mockBuilder, times(1)).build();
    }
  }
}
