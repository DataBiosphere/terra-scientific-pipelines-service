package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.testutils.BaseContainerTest;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class PipelinesServiceTest extends BaseContainerTest {
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
  void testAllPipelineEnumsExist() {
    // make sure all the pipelines in the enum exist in the table
    for (PipelinesEnum p : PipelinesEnum.values()) {
      assertTrue(pipelinesRepository.existsByPipelineId(p.getValue()));
    }
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
