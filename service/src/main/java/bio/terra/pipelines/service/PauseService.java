package bio.terra.pipelines.service;

import bio.terra.pipelines.db.entities.ServiceState;
import bio.terra.pipelines.db.repositories.ServiceStateRepository;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Service;

/**
 * Single source of truth for whether Teaspoons is paused. State lives in the single-row
 * service_state table so it can be flipped without a deploy.
 */
@Service
public class PauseService {

  private static final Logger logger = LoggerFactory.getLogger(PauseService.class);
  private final ServiceStateRepository serviceStateRepository;

  public PauseService(ServiceStateRepository serviceStateRepository) {
    this.serviceStateRepository = serviceStateRepository;
  }

  public boolean isPaused() {
    return getServiceState().isPaused();
  }

  public void pauseService() {
    setPaused(true);
  }

  public void resumeService() {
    setPaused(false);
  }

  private void setPaused(boolean paused) {
    ServiceState state = getServiceState();
    state.setPaused(paused);
    serviceStateRepository.save(state);
    logger.info("Teaspoons service paused state set to {}", paused);
  }

  private ServiceState getServiceState() {
    return serviceStateRepository
        .findById(ServiceState.SINGLETON_ID)
        .orElseThrow(() -> new IllegalStateException("service_state row is missing"));
  }
}
