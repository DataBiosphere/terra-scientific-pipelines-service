package bio.terra.pipelines.common.utils;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.service.ImputationService;
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
  private final ImputationService imputationService;
  private final SamService samService;
  private final LeonardoService leonardoService;
  private final WdsService wdsService;
  private final CbasService cbasService;
  private final ImputationConfiguration imputationConfiguration;

  @Lazy
  @Autowired
  public FlightBeanBag(
      PipelinesService pipelinesService,
      PipelineRunsService pipelineRunsService,
      ImputationService imputationService,
      SamService samService,
      LeonardoService leonardoService,
      WdsService wdsService,
      CbasService cbasService,
      ImputationConfiguration imputationConfiguration) {
    this.pipelinesService = pipelinesService;
    this.pipelineRunsService = pipelineRunsService;
    this.imputationService = imputationService;
    this.samService = samService;
    this.leonardoService = leonardoService;
    this.wdsService = wdsService;
    this.cbasService = cbasService;
    this.imputationConfiguration = imputationConfiguration;
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

  public ImputationService getImputationService() {
    return imputationService;
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
}
