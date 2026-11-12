package bio.terra.pipelines.app.configuration.internal;

import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.retry.RetryLoggingListener;
import bio.terra.rawls.client.ApiException;
import jakarta.ws.rs.ProcessingException;
import java.net.SocketTimeoutException;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.retry.RetryListener;
import org.springframework.retry.annotation.EnableRetry;
import org.springframework.retry.backoff.ExponentialBackOffPolicy;
import org.springframework.retry.backoff.FixedBackOffPolicy;
import org.springframework.retry.policy.ExceptionClassifierRetryPolicy;
import org.springframework.retry.policy.NeverRetryPolicy;
import org.springframework.retry.policy.SimpleRetryPolicy;
import org.springframework.retry.support.RetryTemplate;

/** bean used to retry requests made by the system */
@EnableRetry
@Configuration
public class RetryConfiguration {

  @Bean(name = "listenerResetRetryTemplate")
  public RetryTemplate listenerResetRetryTemplate() {
    RetryTemplate retryTemplate = new RetryTemplate();

    // Fixed delay of 1 second between retries
    FixedBackOffPolicy fixedBackOffPolicy = new FixedBackOffPolicy();
    fixedBackOffPolicy.setBackOffPeriod(1000L);

    // Inner retry (assuming the classifier hits): up to 3 times
    SimpleRetryPolicy srp = new SimpleRetryPolicy();
    srp.setMaxAttempts(3);

    ExceptionClassifierRetryPolicy ecrp = new ExceptionClassifierRetryPolicy();
    ecrp.setExceptionClassifier(
        exception -> {
          if (exception instanceof ProcessingException
              || exception instanceof SocketTimeoutException) {
            return srp;
          } else {
            return new NeverRetryPolicy();
          }
        });

    retryTemplate.setBackOffPolicy(fixedBackOffPolicy);
    retryTemplate.setRetryPolicy(ecrp);
    retryTemplate.setThrowLastExceptionOnExhausted(true);
    retryTemplate.setListeners(new RetryListener[] {new RetryLoggingListener()});

    return retryTemplate;
  }

  /**
   * Retry template for calls to Rawls. In addition to connection-level failures, this retries Rawls
   * 5xx responses (which are often transient) with an exponential backoff for roughly 2 minutes in
   * total: waits of 5s, 10s, 20s, 30s, 30s and 30s between 7 attempts.
   */
  @Bean(name = "rawlsRetryTemplate")
  public RetryTemplate rawlsRetryTemplate() {
    RetryTemplate retryTemplate = new RetryTemplate();

    ExponentialBackOffPolicy backOffPolicy = new ExponentialBackOffPolicy();
    backOffPolicy.setInitialInterval(5000L); // 5 seconds
    backOffPolicy.setMultiplier(2.0); // double the wait time after each retry
    backOffPolicy.setMaxInterval(30000L); // 30 seconds

    // connection-level failures: up to 3 attempts
    SimpleRetryPolicy connectionFailurePolicy = new SimpleRetryPolicy(3);
    // Rawls 5xx responses: up to 7 attempts
    SimpleRetryPolicy serverErrorPolicy = new SimpleRetryPolicy(7);

    ExceptionClassifierRetryPolicy ecrp = new ExceptionClassifierRetryPolicy();
    ecrp.setExceptionClassifier(
        exception -> {
          if (exception instanceof ProcessingException
              || exception instanceof SocketTimeoutException) {
            return connectionFailurePolicy;
          } else if (isRawlsServerError(exception)) {
            return serverErrorPolicy;
          } else {
            return new NeverRetryPolicy();
          }
        });

    retryTemplate.setBackOffPolicy(backOffPolicy);
    retryTemplate.setRetryPolicy(ecrp);
    retryTemplate.setThrowLastExceptionOnExhausted(true);
    retryTemplate.setListeners(new RetryListener[] {new RetryLoggingListener()});

    return retryTemplate;
  }

  /** true if the exception wraps a Rawls ApiException with a 5xx status code */
  private static boolean isRawlsServerError(Throwable exception) {
    return exception instanceof RawlsServiceApiException
        && exception.getCause() instanceof ApiException apiException
        && apiException.getCode() >= 500
        && apiException.getCode() < 600;
  }
}
