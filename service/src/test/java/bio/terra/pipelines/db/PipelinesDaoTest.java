package bio.terra.pipelines.db;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.service.model.Pipeline;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class PipelinesDaoTest extends BaseDaoTest {
  @Autowired PipelinesDao pipelinesDao;

  int nTotalPipelines = 2;

  @Test
  void testGetPipelines() {
    var retrievedPipelines = pipelinesDao.getPipelines();

    assertEquals(retrievedPipelines.size(), nTotalPipelines);

    for (Pipeline pipeline : retrievedPipelines) {
      assertNotNull(pipeline.getPipelineId());
      assertNotNull(pipeline.getDisplayName());
      assertNotNull(pipeline.getDescription());
    }
  }
}
