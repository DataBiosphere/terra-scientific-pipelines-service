package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;

import bio.terra.common.stairway.StairwayComponent;
import bio.terra.pipelines.common.utils.MdcHook;
import bio.terra.pipelines.dependencies.stairway.StairwayJobBuilder;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.testutils.BaseContainerTest;
import java.util.Map;
import org.junit.jupiter.api.BeforeEach;
import org.mockito.InjectMocks;
import org.mockito.Mock;

class JobsServiceMockTest extends BaseContainerTest {

  @Mock private StairwayJobService stairwayJobService;
  @Mock private StairwayComponent stairwayComponent;
  @Mock private MdcHook mdcHook;
  @Mock private StairwayJobBuilder stairwayJobBuilder;

  @InjectMocks private JobsService jobsService;

  // parameters used repeatedly by various tests, and things we'll want mocks to respond to
  // universally
  private final String testUserId = "testUser";
  private final String testGoodPipelineId = "testGoodPipeline";
  private final String testPipelineVersion = "testPipelineVersion";

  private final Object testPipelineInputs = Map.of("first_key", "first_value");
  private final String testUUIDString = "deadbeef-dead-beef-aaaa-beefdeadbeef";

  @BeforeEach
  void initMocks() {
    // stairway submit method returns a good flightId
    when(stairwayJobService.newJob())
        .thenReturn(new StairwayJobBuilder(stairwayJobService, stairwayComponent, mdcHook));
    when(stairwayJobBuilder.jobId(any()))
        .thenReturn(
            new StairwayJobBuilder(stairwayJobService, stairwayComponent, mdcHook)
                .jobId(testUUIDString));
    when(stairwayJobBuilder.submit()).thenReturn(testUUIDString);
  }

  // this is not working and I need help
  //  @Test
  //  void testCreateJob_success() {
  //    String writtenUUID =
  //        jobsService.createJob(
  //            testUserId, testGoodPipelineId, testPipelineVersion, testPipelineInputs);
  //    assertEquals(testUUIDString, writtenUUID);
  //  }
}
