package bio.terra.pipelines.common.utils;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.notifications.NotificationService;
import bio.terra.pipelines.service.DownloadCallCounterService;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.service.QuotasService;
import lombok.Getter;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.Lazy;
import org.springframework.stereotype.Component;

/**
 * The purpose of FlightBeanBag is to provide a clean interface for flights to get access to
 * singleton Spring components. This avoids the use of dynamic bean lookups in flights and casting
 * the lookup result. Instead, flights make calls to accessors in this class. Spring will wire up
 * the underlying methods once at startup avoiding the bean lookup. The objects will be properly
 * types without casting.
 */
@Component
@Getter
public class FlightBeanBag {
  private final PipelinesService pipelinesService;
  private final PipelineRunsService pipelineRunsService;
  private final PipelineInputsOutputsService pipelineInputsOutputsService;
  private final SamService samService;
  private final RawlsService rawlsService;
  private final QuotasService quotasService;
  private final NotificationService notificationService;
  private final DownloadCallCounterService downloadCallCounterService;
  private final PipelineConfigurations pipelineConfigurations;

  @Lazy
  @Autowired
  public FlightBeanBag(
      PipelinesService pipelinesService,
      PipelineRunsService pipelineRunsService,
      PipelineInputsOutputsService pipelineInputsOutputsService,
      SamService samService,
      RawlsService rawlsService,
      QuotasService quotasService,
      NotificationService notificationService,
      DownloadCallCounterService downloadCallCounterService,
      PipelineConfigurations pipelineConfigurations) {
    this.pipelinesService = pipelinesService;
    this.pipelineRunsService = pipelineRunsService;
    this.pipelineInputsOutputsService = pipelineInputsOutputsService;
    this.samService = samService;
    this.rawlsService = rawlsService;
    this.quotasService = quotasService;
    this.notificationService = notificationService;
    this.downloadCallCounterService = downloadCallCounterService;
    this.pipelineConfigurations = pipelineConfigurations;
  }

  public static FlightBeanBag getFromObject(Object object) {
    return (FlightBeanBag) object;
  }
}
