package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.argThat;
import static org.mockito.Mockito.*;

import bio.terra.pipelines.db.entities.Job;
import bio.terra.pipelines.db.repositories.JobsRepository;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.testutils.BaseContainerTest;
import java.util.Map;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;

class JobsServiceMockTest extends BaseContainerTest {

  @InjectMocks private JobsService jobsService;
  @Mock private JobsRepository jobsRepository;
  @Mock private PipelineInputsRepository pipelineInputsRepository;
  private JobsService jobServiceSpy;

  // parameters used repeatedly by various tests, and things we'll want to mocks to respond to
  // universally
  private final String testUserId = "testUser";
  private final String testGoodPipelineId = "testGoodPipeline";
  private final String testPipelineVersion = "testPipelineVersion";

  private final Object testPipelineInputs = Map.of("first_key", "first_value");

  // We'll need these to configure the dao to return selectively good or bad values
  private final UUID testGoodUUID = UUID.randomUUID();
  private final Job testGoodJob = new Job(testGoodUUID, null, null, null, null, null, null);
  private final UUID testDuplicateUUID = UUID.randomUUID();

  @BeforeEach
  void initMocks() {
    // dao returns null on job containing duplicate id and returns good uuid on job containing good
    // uuid
    jobServiceSpy = spy(jobsService);
    doReturn(null)
        .when(jobServiceSpy)
        .writeJobToDbRetryDuplicateException(
            argThat((Job j) -> j.getJobId().equals(testDuplicateUUID)));
    // doReturn is the necessary syntax after an exception-stubbed method.
    // See
    // https://javadoc.io/doc/org.mockito/mockito-core/latest/org/mockito/Mockito.html#doReturn(java.lang.Object)
    doReturn(testGoodJob)
        .when(jobServiceSpy)
        .writeJobToDbRetryDuplicateException(argThat((Job j) -> j.getJobId().equals(testGoodUUID)));
  }

  // JobsService.createJob has 3 pieces of business logic to check.
  // Case 1: It creates a job object.  Can test this by ensuring that it calls the dao with a valid
  // Job object and  doesn't fail if the dao doesn't.  Beyond not returning a null, verify that it
  // returns the same UUID
  // Case 2: If a duplicate key error occurs consistently while writing the job object, it should
  // return null
  // Case 3: If a duplicate key error occurs once while writing the job object, it should retry and
  // return the successfully written UUID
  @Test
  void testCreateJob_successfulWriteUUIDsMatch() {
    // override a bit of our bean with a spy here, which leaves the rest untouched
    doReturn(testGoodUUID).when(jobServiceSpy).createJobId();

    UUID writtenUUID =
        jobServiceSpy.createJob(
            testUserId, testGoodPipelineId, testPipelineVersion, testPipelineInputs);
    assertEquals(writtenUUID, testGoodUUID);
  }

  @Test
  void testCreateJob_unsuccessfulWriteDaoReturnsNull() {
    // override a bit of our bean with a spy here, which leaves the rest untouched
    doReturn(testDuplicateUUID).when(jobServiceSpy).createJobId();
    UUID returnedUUID =
        jobServiceSpy.createJob(
            testUserId, testGoodPipelineId, testPipelineVersion, testPipelineInputs);

    assertNull(returnedUUID);
  }

  @Test
  void testCreateJob_unsuccessfulWriteDaoReturnsNullThenSucceeds() {
    // override a bit of our bean with a spy here, which leaves the rest untouched
    doReturn(testDuplicateUUID, testGoodUUID).when(jobServiceSpy).createJobId();
    UUID returnedUUID =
        jobServiceSpy.createJob(
            testUserId, testGoodPipelineId, testPipelineVersion, testPipelineInputs);

    assertEquals(returnedUUID, testGoodUUID);
  }
}
