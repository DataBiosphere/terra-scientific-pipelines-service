package bio.terra.pipelines.retry;

import bio.terra.pipelines.app.configuration.internal.RetryConfiguration;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.retry.RetryCallback;
import org.springframework.retry.RetryContext;
import org.springframework.retry.RetryListener;
import org.springframework.stereotype.Component;

/**
 * Listener used to log messages on requests that are retried through {@link
 * RetryConfiguration#listenerResetRetryTemplate}
 */
@Component("retryLoggingListener")
public class RetryLoggingListener implements RetryListener {

  private final Logger logger = LoggerFactory.getLogger(getClass());

  @Override
  public <T, E extends Throwable> void close(
      RetryContext context, RetryCallback<T, E> callback, Throwable throwable) {
    if (context.getRetryCount() > 1) {
      logger.debug("Retryable method closing ({}th retry)", context.getRetryCount());
    }
  }

  @Override
  public <T, E extends Throwable> boolean open(RetryContext context, RetryCallback<T, E> callback) {
    if (context.getRetryCount() > 1) {
      logger.debug("Retryable method opening (retry count: {})", context.getRetryCount());
    }
    return true;
  }

  @Override
  public <T, E extends Throwable> void onError(
      RetryContext context, RetryCallback<T, E> callback, Throwable throwable) {
    logger.warn(
        "Retryable method threw exception (retry count: {})", context.getRetryCount(), throwable);
  }
}
