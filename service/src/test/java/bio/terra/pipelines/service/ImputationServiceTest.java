package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.Map;
import java.util.UUID;
import org.springframework.beans.factory.annotation.Autowired;

class ImputationServiceTest extends BaseEmbeddedDbTest {

  @Autowired ImputationService imputationService;

  private final String testUserId = TestUtils.TEST_USER_ID_1;

  private final String testStatus = "RUNNING";
  private final Long testPipelineId = TestUtils.TEST_PIPELINE_ID_1;
  private final Map<String, Object> testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testJobId = TestUtils.TEST_NEW_UUID;
}
