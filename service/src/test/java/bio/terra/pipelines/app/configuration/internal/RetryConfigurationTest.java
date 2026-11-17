package bio.terra.pipelines.app.configuration.internal;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;

import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.rawls.client.ApiException;
import jakarta.ws.rs.ProcessingException;
import java.net.SocketTimeoutException;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.retry.backoff.FixedBackOffPolicy;
import org.springframework.retry.support.RetryTemplate;

/**
 * Unit tests for which exceptions the rawlsRetryTemplate retries, and how many attempts it makes.
 */
class RetryConfigurationTest {

  private static final int SERVER_ERROR_MAX_ATTEMPTS = 7;
  private static final int CONNECTION_FAILURE_MAX_ATTEMPTS = 3;

  private RetryTemplate template;

  @BeforeEach
  void setup() {
    template = new RetryConfiguration().rawlsRetryTemplate();
    // use a tiny backoff so the tests don't actually wait
    FixedBackOffPolicy smallBackoff = new FixedBackOffPolicy();
    smallBackoff.setBackOffPeriod(1L);
    template.setBackOffPolicy(smallBackoff);
  }

  /**
   * Run a call that always throws the given exception through the template and return how many
   * times it was attempted. Also checks that the original exception is the one that is rethrown.
   */
  private int attemptsFor(Exception exceptionToThrow) {
    AtomicInteger attempts = new AtomicInteger();
    Exception thrown =
        assertThrows(
            Exception.class,
            () ->
                template.<Object, Exception>execute(
                    context -> {
                      attempts.incrementAndGet();
                      throw exceptionToThrow;
                    }));
    assertSame(exceptionToThrow, thrown);
    return attempts.get();
  }

  private static RawlsServiceApiException rawlsApiException(int statusCode) {
    return new RawlsServiceApiException(new ApiException(statusCode, "status " + statusCode));
  }

  @Test
  void serverErrorsAreRetriedUpToMaxAttempts() {
    for (int statusCode : List.of(500, 502, 503, 504, 599)) {
      assertEquals(
          SERVER_ERROR_MAX_ATTEMPTS,
          attemptsFor(rawlsApiException(statusCode)),
          "status code " + statusCode);
    }
  }

  @Test
  void nonServerErrorStatusCodesAreNotRetried() {
    for (int statusCode : List.of(400, 404, 409, 429, 499, 600)) {
      assertEquals(1, attemptsFor(rawlsApiException(statusCode)), "status code " + statusCode);
    }
  }

  @Test
  void connectionFailuresAreRetriedUpToMaxAttempts() {
    assertEquals(
        CONNECTION_FAILURE_MAX_ATTEMPTS, attemptsFor(new ProcessingException("connection reset")));
    assertEquals(
        CONNECTION_FAILURE_MAX_ATTEMPTS, attemptsFor(new SocketTimeoutException("timed out")));
  }

  @Test
  void otherExceptionsAreNotRetried() {
    assertEquals(1, attemptsFor(new IllegalStateException("unexpected")));
    // a RawlsServiceApiException that is not caused by an ApiException is not a Rawls server error
    assertEquals(1, attemptsFor(new RawlsServiceApiException("no api exception cause")));
  }

  @Test
  void serverErrorRetriesEventuallySucceed() throws Exception {
    AtomicInteger attempts = new AtomicInteger();
    String result =
        template.execute(
            context -> {
              if (attempts.incrementAndGet() < 4) {
                throw rawlsApiException(500);
              }
              return "success";
            });
    assertEquals("success", result);
    assertEquals(4, attempts.get());
  }
}
