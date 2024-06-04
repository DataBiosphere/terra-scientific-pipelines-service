package bio.terra.pipelines.common.utils;

import bio.terra.pipelines.app.configuration.external.CbasConfiguration;
import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.PipelinesService;
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
public class FlightBeanBag {
  private final PipelinesService pipelinesService;
  private final PipelineRunsService pipelineRunsService;
  private final SamService samService;
  private final LeonardoService leonardoService;
  private final WdsService wdsService;
  private final CbasService cbasService;
  private final ImputationConfiguration imputationConfiguration;
  private final CbasConfiguration cbasConfiguration;

  @Lazy
  @Autowired
  public FlightBeanBag(
      PipelinesService pipelinesService,
      PipelineRunsService pipelineRunsService,
      SamService samService,
      LeonardoService leonardoService,
      WdsService wdsService,
      CbasService cbasService,
      ImputationConfiguration imputationConfiguration,
      CbasConfiguration cbasConfiguration) {
    this.pipelinesService = pipelinesService;
    this.pipelineRunsService = pipelineRunsService;
    this.samService = samService;
    this.leonardoService = leonardoService;
    this.wdsService = wdsService;
    this.cbasService = cbasService;
    this.imputationConfiguration = imputationConfiguration;
    this.cbasConfiguration = cbasConfiguration;
  }

  public static FlightBeanBag getFromObject(Object object) {
    return (FlightBeanBag) object;
  }

  public PipelinesService getPipelinesService() {
    return pipelinesService;
  }

  public PipelineRunsService getPipelineRunsService() {
    return pipelineRunsService;
  }

  public SamService getSamService() {
    return samService;
  }

  public LeonardoService getLeonardoService() {
    return leonardoService;
  }

  public WdsService getWdsService() {
    return wdsService;
  }

  public CbasService getCbasService() {
    return cbasService;
  }

  public ImputationConfiguration getImputationConfiguration() {
    return imputationConfiguration;
  }

  public CbasConfiguration getCbasConfiguration() {
    return cbasConfiguration;
  }
}
