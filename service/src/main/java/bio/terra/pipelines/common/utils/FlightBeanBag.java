package bio.terra.pipelines.common.utils;

import bio.terra.pipelines.service.ImputationService;
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
  private final ImputationService imputationService;

  @Lazy
  @Autowired
  public FlightBeanBag(ImputationService imputationService) {
    this.imputationService = imputationService;
  }

  public static FlightBeanBag getFromObject(Object object) {
    return (FlightBeanBag) object;
  }

  public ImputationService getImputationService() {
    return imputationService;
  }
}
