package bio.terra.pipelines.common.utils;

import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.dependencies.sam.SamService;
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
  private final PipelinesRepository pipelinesRepository;
  private final SamService samService;
  private final PipelinesService pipelinesService;

  @Lazy
  @Autowired
  public FlightBeanBag(
      PipelinesRepository pipelinesRepository,
      SamService samService,
      PipelinesService pipelinesService) {
    this.pipelinesRepository = pipelinesRepository;
    this.samService = samService;
    this.pipelinesService = pipelinesService;
  }

  public static FlightBeanBag getFromObject(Object object) {
    return (FlightBeanBag) object;
  }

  public PipelinesRepository getPipelinesRepository() {
    return pipelinesRepository;
  }

  public SamService getSamService() {
    return samService;
  }

  public PipelinesService getPipelinesService() {
    return pipelinesService;
  }
}
