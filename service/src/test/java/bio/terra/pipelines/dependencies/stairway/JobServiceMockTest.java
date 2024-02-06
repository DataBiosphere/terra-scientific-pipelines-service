package bio.terra.pipelines.dependencies.stairway;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.common.stairway.StairwayComponent;
import bio.terra.pipelines.app.controller.JobApiUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.exception.*;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.MakeFlightException;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.Optional;
import java.util.UUID;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.ArgumentMatchers;
import org.mockito.InjectMocks;
import org.mockito.Mock;

class JobServiceMockTest extends BaseEmbeddedDbTest {

  @InjectMocks JobService jobService;
  @Mock private Stairway mockStairway;
  @Mock private StairwayComponent mockStairwayComponent;

  @BeforeEach
  void setup() {
    when(mockStairwayComponent.get()).thenReturn(mockStairway);
  }

  /**
   * Reset the {@link JobService} {@link FlightDebugInfo} after each test so that future submissions
   * aren't affected.
   */
  @AfterEach
  void clearFlightDebugInfo() {
    jobService.setFlightDebugInfoForTest(null);
  }

  @Test
  void submit_StairwayException() throws InterruptedException {
    // a MakeFlightException( is an instance of StairwayException
    doThrow(new MakeFlightException("test exception"))
        .when(mockStairway)
        .submitWithDebugInfo(any(), any(), any(), ArgumentMatchers.eq(false), any());

    UUID jobId = TestUtils.TEST_NEW_UUID;
    assertThrows(InternalStairwayException.class, () -> jobService.submit(null, null, jobId));
  }

  @Test
  void retrieveJobResult_successWithResultClass() throws InterruptedException {
    FlightMap inputParams = new FlightMap();
    FlightMap flightMap = new FlightMap();
    String expectedResponse = "foo";
    flightMap.put(JobMapKeys.RESPONSE.getKeyName(), expectedResponse);
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState successFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.SUCCESS, flightId, inputParams, flightMap);

    when(mockStairway.getFlightState(any())).thenReturn(successFlightState);

    JobService.JobResultOrException<String> resultOrException =
        jobService.retrieveJobResult(flightId, String.class);
    assertEquals(expectedResponse, resultOrException.getResult());
  }

  @Test
  void retrieveJobResult_successWithResultTypeRef() throws InterruptedException {
    FlightMap inputParams = new FlightMap();
    FlightMap flightMap = new FlightMap();
    String expectedResponse = "foo";
    flightMap.put(JobMapKeys.RESPONSE.getKeyName(), expectedResponse);
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState successFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.SUCCESS, flightId, inputParams, flightMap);

    when(mockStairway.getFlightState(any())).thenReturn(successFlightState);

    JobService.JobResultOrException<String> resultOrException =
        jobService.retrieveJobResult(flightId, null, new TypeReference<>() {});
    assertEquals(expectedResponse, resultOrException.getResult());
  }

  @Test
  void retrieveJobResult_noResultClassOrTypeThrows() throws InterruptedException {
    FlightMap inputParams = new FlightMap();
    FlightMap flightMap = new FlightMap();
    flightMap.put(JobMapKeys.RESPONSE.getKeyName(), null);
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState successFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.SUCCESS, flightId, inputParams, flightMap);

    when(mockStairway.getFlightState(any())).thenReturn(successFlightState);

    assertThrows(
        InvalidResultStateException.class,
        () -> jobService.retrieveJobResult(flightId, null, null));
  }

  @Test
  void retrieveJobResult_fatal() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState fatalFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.FATAL, flightId);
    fatalFlightState.setException(new RuntimeException("test exception"));

    when(mockStairway.getFlightState(any())).thenReturn(fatalFlightState);

    JobService.JobResultOrException result =
        jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class);
    assertEquals(fatalFlightState.getException(), Optional.ofNullable(result.getException()));
  }

  @Test
  void retrieveJobResult_fatalNonRuntime() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState fatalFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.FATAL, flightId);
    // non-runtime exception should be caught and stored as an InternalServerErrorException
    fatalFlightState.setException(new InterruptedException());

    when(mockStairway.getFlightState(any())).thenReturn(fatalFlightState);

    JobService.JobResultOrException result =
        jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class);
    assertInstanceOf(InternalServerErrorException.class, result.getException());
  }

  @Test
  void retrieveJobResult_errorFlightState() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState errorFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.ERROR, flightId);
    errorFlightState.setException(new RuntimeException("test exception"));

    when(mockStairway.getFlightState(any())).thenReturn(errorFlightState);

    JobService.JobResultOrException result =
        jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class);
    assertEquals(errorFlightState.getException(), Optional.ofNullable(result.getException()));
  }

  @Test
  void retrieveJobResult_StairwayException() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    // a MakeFlightException( is an instance of StairwayException
    when(mockStairway.getFlightState(any())).thenThrow(new MakeFlightException("test exception"));

    assertThrows(
        InternalStairwayException.class,
        () -> jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class));
  }

  @Test
  void retrieveJobResult_running() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState runningFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.RUNNING, flightId);

    when(mockStairway.getFlightState(any())).thenReturn(runningFlightState);

    assertThrows(
        JobNotCompleteException.class,
        () -> jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class, null));
  }

  @Test
  void retrieveJobResult_interrupted() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    when(mockStairway.getFlightState(any())).thenThrow(new InterruptedException("test exception"));

    assertThrows(
        InternalStairwayException.class,
        () -> jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class, null));
  }

  @Test
  void retrieveJob_stairwayException() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    String userId = "testUserId";
    // a MakeFlightException( is an instance of StairwayException
    when(mockStairway.getFlightState(any())).thenThrow(new MakeFlightException("test exception"));

    assertThrows(
        InternalStairwayException.class, () -> jobService.retrieveJob(flightId, userId, null));
  }

  @Test
  void retrieveJob_interruptedException() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    String userId = "testUserId";

    when(mockStairway.getFlightState(any())).thenThrow(new InterruptedException());

    // InterruptedException should be caught and re-thrown as an InternalStairwayException
    assertThrows(
        InternalStairwayException.class, () -> jobService.retrieveJob(flightId, userId, null));
  }

  @Test
  void testRetrieveAsyncJobResultRunning() throws InterruptedException {
    UUID jobId = TestUtils.TEST_NEW_UUID;
    FlightMap inputParameters = StairwayTestUtils.CREATE_JOB_INPUT_PARAMS;
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.RUNNING, jobId, inputParameters, new FlightMap());

    when(mockStairway.getFlightState(any())).thenReturn(flightState);

    JobApiUtils.AsyncJobResult<String> result =
        jobService.retrieveAsyncJobResult(
            jobId, TestUtils.TEST_USER_ID_1, PipelinesEnum.IMPUTATION_MINIMAC4, String.class, null);

    assertEquals(jobId.toString(), result.getJobReport().getId());
    assertEquals(202, result.getJobReport().getStatusCode());
    assertNull(result.getResult());
    assertNull(result.getApiErrorReport());
  }

  @Test
  void testRetrieveAsyncJobResultSucceeded() throws InterruptedException {
    UUID jobId = TestUtils.TEST_NEW_UUID;
    FlightMap inputParameters = StairwayTestUtils.CREATE_JOB_INPUT_PARAMS;
    FlightMap workingMap = new FlightMap();
    String testResponse = "test response";
    workingMap.put(JobMapKeys.RESPONSE.getKeyName(), testResponse);
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.SUCCESS, jobId, inputParameters, workingMap);

    when(mockStairway.getFlightState(any())).thenReturn(flightState);

    JobApiUtils.AsyncJobResult<String> result =
        jobService.retrieveAsyncJobResult(
            jobId, TestUtils.TEST_USER_ID_1, PipelinesEnum.IMPUTATION_MINIMAC4, String.class, null);

    assertEquals(jobId.toString(), result.getJobReport().getId());
    assertEquals(200, result.getJobReport().getStatusCode());
    assertEquals(testResponse, result.getResult());
    assertNull(result.getApiErrorReport());
  }

  @Test
  void testRetrieveAsyncJobResultFailed() throws InterruptedException {
    UUID jobId = TestUtils.TEST_NEW_UUID;
    FlightMap inputParameters = StairwayTestUtils.CREATE_JOB_INPUT_PARAMS;
    // even on a fatal failure the response might have been written to the working map
    FlightMap workingMap = new FlightMap();
    String testResponse = "test response";
    workingMap.put(JobMapKeys.RESPONSE.getKeyName(), testResponse);
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.ERROR, jobId, inputParameters, workingMap);
    String testErrorMsg = "test exception";
    flightState.setException(new RuntimeException(testErrorMsg));

    when(mockStairway.getFlightState(any())).thenReturn(flightState);

    JobApiUtils.AsyncJobResult<String> result =
        jobService.retrieveAsyncJobResult(
            jobId, TestUtils.TEST_USER_ID_1, PipelinesEnum.IMPUTATION_MINIMAC4, String.class, null);

    assertEquals(jobId.toString(), result.getJobReport().getId());
    assertEquals(500, result.getJobReport().getStatusCode());
    assertNull(result.getResult());
    assertEquals(testErrorMsg, result.getApiErrorReport().getMessage());
  }

  @Test
  void retrieveAsyncJobResultStairwayException() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    // a MakeFlightException is an instance of StairwayException
    when(mockStairway.getFlightState(any())).thenThrow(new MakeFlightException("test exception"));

    assertThrows(
        InternalStairwayException.class,
        () ->
            jobService.retrieveAsyncJobResult(
                flightId,
                TestUtils.TEST_USER_ID_1,
                TestUtils.TEST_PIPELINE_1_ENUM,
                String.class,
                null));
  }

  @Test
  void enumerateJobs_stairwayException() throws InterruptedException {
    String userId = "testUserId";
    // a MakeFlightException( is an instance of StairwayException
    when(mockStairway.getFlights(any(), any(), any()))
        .thenThrow(new MakeFlightException("test exception"));

    assertThrows(
        InternalStairwayException.class, () -> jobService.enumerateJobs(userId, 10, null, null));
  }

  @Test
  void enumerateJobs_interruptedException() throws InterruptedException {
    String userId = "testUserId";

    when(mockStairway.getFlights(any(), any(), any())).thenThrow(new InterruptedException());

    // InterruptedException should be caught and re-thrown as an InternalStairwayException
    assertThrows(
        InternalStairwayException.class, () -> jobService.enumerateJobs(userId, 10, null, null));
  }
}
