package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.db.entities.ServiceState;
import bio.terra.pipelines.db.repositories.ServiceStateRepository;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import java.time.Instant;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class PauseServiceTest extends BaseEmbeddedDbTest {

  @Autowired PauseService pauseService;
  @Autowired ServiceStateRepository serviceStateRepository;

  @Test
  void seededRowIsUnpaused() {
    assertEquals(1, serviceStateRepository.count());
    assertFalse(pauseService.isPaused());
  }

  @Test
  void pauseAndResumeService() {
    Instant before = serviceStateRepository.findById(ServiceState.SINGLETON_ID).get().getUpdated();

    pauseService.pauseService();
    assertTrue(pauseService.isPaused());
    Instant afterPause =
        serviceStateRepository.findById(ServiceState.SINGLETON_ID).get().getUpdated();
    assertTrue(afterPause.isAfter(before));

    pauseService.resumeService();
    assertFalse(pauseService.isPaused());
    assertEquals(1, serviceStateRepository.count());
  }
}
