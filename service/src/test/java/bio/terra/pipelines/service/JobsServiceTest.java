package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.argThat;
import static org.mockito.Mockito.*;

import bio.terra.pipelines.db.JobsDao;
import bio.terra.pipelines.service.model.Job;
import bio.terra.pipelines.testutils.BaseUnitTest;
import java.util.Map;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;

class JobsServiceTest extends BaseUnitTest {

  @Autowired private JobsService jobsService;
  @MockBean private JobsDao jobsDao;

  // parameters used repeatedly by various tests, and things we'll want to mocks to respond to
  // universally
  private final String testUserId = "testUser";
  private final String testGoodPipelineId = "testGoodPipeline";
  private final String testPipelineVersion = "testPipelineVersion";

  private final Object testPipelineInputs = Map.of();

  // We'll need these to configure the dao to return selectively good or bad values
  private final UUID testGoodUUID = UUID.randomUUID();
  private final UUID testDuplicateUUID = UUID.randomUUID();

  @BeforeEach
  void initMocks() {
    // dao returns null on job containing duplicate id and returns good uuid on job containing good
    // uuid
    when(jobsDao.createJob(argThat((Job j) -> j.getJobId() == testDuplicateUUID))).thenReturn(null);
    // doReturn is the necessary syntax after an exception-stubbed method.
    // See:
    // https://javadoc.io/doc/org.mockito/mockito-core/latest/org/mockito/Mockito.html#doReturn(java.lang.Object)
    doReturn(testGoodUUID)
        .when(jobsDao)
        .createJob(argThat((Job j) -> j.getJobId() == testGoodUUID));
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
    JobsService jobServiceSpy = spy(jobsService);
    doReturn(testGoodUUID).when(jobServiceSpy).createJobId();

    UUID writtenUUID =
        jobServiceSpy.createJob(
            testUserId, testGoodPipelineId, testPipelineVersion, testPipelineInputs);
    assertEquals(writtenUUID, testGoodUUID);
  }

  @Test
  void testCreateJob_unsuccessfulWriteDaoReturnsNull() {
    // override a bit of our bean with a spy here, which leaves the rest untouched
    JobsService jobServiceSpy = spy(jobsService);
    doReturn(testDuplicateUUID).when(jobServiceSpy).createJobId();
    UUID returnedUUID =
        jobServiceSpy.createJob(
            testUserId, testGoodPipelineId, testPipelineVersion, testPipelineInputs);

    assertNull(returnedUUID);
  }

  @Test
  void testCreateJob_unsuccessfulWriteDaoReturnsNullThenSucceeds() {
    // override a bit of our bean with a spy here, which leaves the rest untouched
    JobsService jobServiceSpy = spy(jobsService);
    doReturn(testDuplicateUUID, testGoodUUID).when(jobServiceSpy).createJobId();
    UUID returnedUUID =
        jobServiceSpy.createJob(
            testUserId, testGoodPipelineId, testPipelineVersion, testPipelineInputs);

    assertEquals(returnedUUID, testGoodUUID);
  }
}
