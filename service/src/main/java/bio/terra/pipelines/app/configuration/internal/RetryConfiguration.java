package bio.terra.pipelines.app.configuration.internal;

import bio.terra.pipelines.retry.RetryLoggingListener;
import jakarta.ws.rs.ProcessingException;
import java.net.SocketTimeoutException;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.retry.RetryListener;
import org.springframework.retry.annotation.EnableRetry;
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
}
