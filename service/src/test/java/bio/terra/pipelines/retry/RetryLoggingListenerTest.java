package bio.terra.pipelines.retry;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;
import ch.qos.logback.classic.spi.ILoggingEvent;
import ch.qos.logback.core.read.ListAppender;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.slf4j.LoggerFactory;
import org.springframework.retry.RetryListener;
import org.springframework.retry.backoff.FixedBackOffPolicy;
import org.springframework.retry.policy.SimpleRetryPolicy;
import org.springframework.retry.support.RetryTemplate;

class RetryLoggingListenerTest {

  private final Logger listenerLogger =
      (Logger) LoggerFactory.getLogger(RetryLoggingListener.class);
  private final ListAppender<ILoggingEvent> logAppender = new ListAppender<>();
  private Level originalLevel;

  private RetryTemplate template;

  @BeforeEach
  void setup() {
    originalLevel = listenerLogger.getLevel();
    listenerLogger.setLevel(Level.DEBUG);
    logAppender.start();
    listenerLogger.addAppender(logAppender);

    template = new RetryTemplate();
    template.setRetryPolicy(new SimpleRetryPolicy(3));
    FixedBackOffPolicy smallBackoff = new FixedBackOffPolicy();
    smallBackoff.setBackOffPeriod(1L);
    template.setBackOffPolicy(smallBackoff);
    template.setThrowLastExceptionOnExhausted(true);
    template.setListeners(new RetryListener[] {new RetryLoggingListener()});
  }

  @AfterEach
  void tearDown() {
    listenerLogger.detachAppender(logAppender);
    logAppender.stop();
    listenerLogger.setLevel(originalLevel);
  }

  private List<ILoggingEvent> warnEvents() {
    return logAppender.list.stream().filter(event -> event.getLevel() == Level.WARN).toList();
  }

  @Test
  void onErrorLogsIncrementingFailedAttemptCount() {
    AtomicInteger attempts = new AtomicInteger();
    assertThrows(
        IllegalStateException.class,
        () ->
            template.<Object, IllegalStateException>execute(
                context -> {
                  attempts.incrementAndGet();
                  throw new IllegalStateException("boom");
                }));
    assertEquals(3, attempts.get());

    List<ILoggingEvent> warnings = warnEvents();
    assertEquals(3, warnings.size());
    for (int i = 0; i < warnings.size(); i++) {
      ILoggingEvent event = warnings.get(i);
      assertTrue(
          event.getFormattedMessage().contains("failed attempt count: " + (i + 1)),
          "unexpected message: " + event.getFormattedMessage());
      // the exception is logged along with the count
      assertNotNull(event.getThrowableProxy());
    }
  }

  @Test
  void nothingLoggedAtWarnWhenFirstAttemptSucceeds() {
    String result = template.execute(context -> "success");
    assertEquals("success", result);
    assertTrue(warnEvents().isEmpty());
  }
}
