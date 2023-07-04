package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.orm.jpa.DataJpaTest;
import org.springframework.test.context.ActiveProfiles;

@DataJpaTest(properties = "spring.main.lazy-initialization=true")
@ActiveProfiles({"test", "human-readable-logging"})
class PipelinesServiceTest {
  @Autowired PipelinesService pipelinesService;
  @Autowired PipelinesRepository pipelinesRepository;

  @Test
  void testGetCorrectNumberOfPipelines() {
    // migrations insert two different pipelines so make sure we find those
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    assertEquals(2, pipelineList.size());

    pipelinesRepository.save(
        new Pipeline("pipelineId", "1.0.0", "pipelineDisplayName", "description"));

    pipelineList = pipelinesService.getPipelines();
    assertEquals(3, pipelineList.size());
  }

  @Test
  void testPipelineExists() {
    // migrations insert an imputation pipeline so that should already exist in the table and the
    // other should not
    assertTrue(pipelinesService.pipelineExists("imputation"));
    assertFalse(pipelinesService.pipelineExists("doesnotexist"));
  }

  @Test
  void testPipelineToString() {
    // test .ToString() method on Pipeline Entity
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    assertEquals(2, pipelineList.size());
    for (Pipeline p : pipelineList) {
      assertEquals(
          String.format(
              "Pipeline[pipelineId=%s, version=%s, displayName=%s, description=%s]",
              p.getPipelineId(), p.getVersion(), p.getDisplayName(), p.getDescription()),
          p.toString());
    }
  }
}
