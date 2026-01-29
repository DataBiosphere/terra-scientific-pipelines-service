package bio.terra.pipelines.service;

import bio.terra.pipelines.app.configuration.internal.StatusCheckConfiguration;
import bio.terra.pipelines.generated.model.ApiSystemStatusSystems;
import com.google.common.annotations.VisibleForTesting;
import jakarta.annotation.PostConstruct;
import java.time.Instant;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class BaseStatusService {
  private static final Logger logger = LoggerFactory.getLogger(BaseStatusService.class);

  /** cached status */
  private final AtomicBoolean cachedStatus;

  /** configuration parameters */
  private final StatusCheckConfiguration configuration;

  /** set of status methods to check */
  private final ConcurrentHashMap<String, Supplier<ApiSystemStatusSystems>> statusCheckMap;

  /** scheduler */
  private final ScheduledExecutorService scheduler;

  /** last time cache was updated */
  private final AtomicReference<Instant> lastStatusUpdate;

  public BaseStatusService(StatusCheckConfiguration configuration) {
    this.configuration = configuration;
    statusCheckMap = new ConcurrentHashMap<>();
    cachedStatus = new AtomicBoolean(false);
    lastStatusUpdate = new AtomicReference<>(Instant.now());
    scheduler = Executors.newScheduledThreadPool(1);
  }

  @PostConstruct
  private void startStatusChecking() {
    if (configuration.enabled()) {
      scheduler.scheduleAtFixedRate(
          this::checkStatus,
          configuration.startupWaitSeconds(),
          configuration.pollingIntervalSeconds(),
          TimeUnit.SECONDS);
    }
  }

  void registerStatusCheck(String name, Supplier<ApiSystemStatusSystems> checkFn) {
    statusCheckMap.put(name, checkFn);
  }

  @VisibleForTesting
  void checkStatus() {
    if (configuration.enabled()) {
      AtomicBoolean newStatus = new AtomicBoolean(true);
      try {
        var systems =
            statusCheckMap.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, e -> e.getValue().get()));
        newStatus.set(systems.values().stream().allMatch(ApiSystemStatusSystems::isOk));
      } catch (Exception e) {
        logger.warn("Status check exception", e);
        newStatus.set(false);
      }
      cachedStatus.set(newStatus.get());
      lastStatusUpdate.set(Instant.now());
    }
  }

  public boolean getCurrentStatus() {
    if (configuration.enabled()) {
      // If staleness time (last update + stale threshold) is before the current time, then
      // we are officially not OK.
      if (lastStatusUpdate
          .get()
          .plusSeconds(configuration.stalenessThresholdSeconds())
          .isBefore(Instant.now())) {

        logger.warn(
            "Status has not been updated since {}. This might mean that the status cronjob has failed, or that requests to downstream services are timing out.",
            lastStatusUpdate);
        cachedStatus.set(false);
      }
      return cachedStatus.get();
    }
    return true;
  }
}
