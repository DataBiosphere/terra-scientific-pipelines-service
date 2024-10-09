package bio.terra.pipelines.dependencies.stairway;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.common.stairway.StairwayComponent;
import bio.terra.pipelines.app.controller.JobApiUtils;
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
  void submitStairwayException() throws InterruptedException {
    UUID jobId = TestUtils.TEST_NEW_UUID;
    // a MakeFlightException is an instance of StairwayException
    doThrow(new MakeFlightException("test exception"))
        .when(mockStairway)
        .submitWithDebugInfo(jobId.toString(), null, null, false, null);

    assertThrows(InternalStairwayException.class, () -> jobService.submit(null, null, jobId));
  }

  @Test
  void retrieveJobResultSuccessWithResultClass() throws InterruptedException {
    FlightMap inputParams = new FlightMap();
    FlightMap flightMap = new FlightMap();
    String expectedResponse = "foo";
    flightMap.put(JobMapKeys.RESPONSE, expectedResponse);
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState successFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.SUCCESS, flightId, inputParams, flightMap);

    when(mockStairway.getFlightState(flightId.toString())).thenReturn(successFlightState);

    JobService.JobResultOrException<String> resultOrException =
        jobService.retrieveJobResult(flightId, String.class);
    assertEquals(expectedResponse, resultOrException.getResult());
  }

  @Test
  void retrieveJobResultSuccessWithResultTypeRef() throws InterruptedException {
    FlightMap inputParams = new FlightMap();
    FlightMap flightMap = new FlightMap();
    String expectedResponse = "foo";
    flightMap.put(JobMapKeys.RESPONSE, expectedResponse);
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState successFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.SUCCESS, flightId, inputParams, flightMap);

    when(mockStairway.getFlightState(flightId.toString())).thenReturn(successFlightState);

    JobService.JobResultOrException<String> resultOrException =
        jobService.retrieveJobResult(flightId, null, new TypeReference<>() {});
    assertEquals(expectedResponse, resultOrException.getResult());
  }

  @Test
  void retrieveJobResultNoResultClassOrTypeThrows() throws InterruptedException {
    FlightMap inputParams = new FlightMap();
    FlightMap flightMap = new FlightMap();
    flightMap.put(JobMapKeys.RESPONSE, null);
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState successFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.SUCCESS, flightId, inputParams, flightMap);

    when(mockStairway.getFlightState(flightId.toString())).thenReturn(successFlightState);

    assertThrows(
        InvalidResultStateException.class,
        () -> jobService.retrieveJobResult(flightId, null, null));
  }

  @Test
  void retrieveJobResultFatal() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState fatalFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.FATAL, flightId);
    fatalFlightState.setException(new RuntimeException("test exception"));

    when(mockStairway.getFlightState(flightId.toString())).thenReturn(fatalFlightState);

    JobService.JobResultOrException result =
        jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class);
    assertEquals(fatalFlightState.getException(), Optional.ofNullable(result.getException()));
  }

  @Test
  void retrieveJobResultFatalNonRuntime() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState fatalFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.FATAL, flightId);
    // non-runtime exception should be caught and stored as an InternalServerErrorException
    fatalFlightState.setException(new InterruptedException());

    when(mockStairway.getFlightState(flightId.toString())).thenReturn(fatalFlightState);

    JobService.JobResultOrException result =
        jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class);
    assertInstanceOf(InternalServerErrorException.class, result.getException());
  }

  @Test
  void retrieveJobResultErrorFlightState() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState errorFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.ERROR, flightId);
    errorFlightState.setException(new RuntimeException("test exception"));

    when(mockStairway.getFlightState(flightId.toString())).thenReturn(errorFlightState);

    JobService.JobResultOrException result =
        jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class);
    assertEquals(errorFlightState.getException(), Optional.ofNullable(result.getException()));
  }

  @Test
  void retrieveJobResultStairwayException() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    // a MakeFlightException is an instance of StairwayException
    when(mockStairway.getFlightState(flightId.toString()))
        .thenThrow(new MakeFlightException("test exception"));

    assertThrows(
        InternalStairwayException.class,
        () -> jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class));
  }

  @Test
  void retrieveJobResultRunning() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    FlightState runningFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.RUNNING, flightId);

    when(mockStairway.getFlightState(flightId.toString())).thenReturn(runningFlightState);

    assertThrows(
        JobNotCompleteException.class,
        () -> jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class, null));
  }

  @Test
  void retrieveJobResultInterrupted() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    when(mockStairway.getFlightState(flightId.toString()))
        .thenThrow(new InterruptedException("test exception"));

    assertThrows(
        InternalStairwayException.class,
        () -> jobService.retrieveJobResult(flightId, JobService.JobResultOrException.class, null));
  }

  @Test
  void retrieveJobStairwayException() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    String userId = "testUserId";
    // a MakeFlightException is an instance of StairwayException
    when(mockStairway.getFlightState(flightId.toString()))
        .thenThrow(new MakeFlightException("test exception"));

    assertThrows(
        InternalStairwayException.class, () -> jobService.retrieveJob(flightId, userId, null));
  }

  @Test
  void retrieveJobInterruptedException() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    String userId = "testUserId";

    when(mockStairway.getFlightState(flightId.toString())).thenThrow(new InterruptedException());

    // InterruptedException should be caught and re-thrown as an InternalStairwayException
    assertThrows(
        InternalStairwayException.class, () -> jobService.retrieveJob(flightId, userId, null));
  }

  @Test
  void enumerateJobsStairwayException() throws InterruptedException {
    String userId = "testUserId";
    // a MakeFlightException is an instance of StairwayException
    when(mockStairway.getFlights((String) eq(null), eq(10), any(FlightFilter.class)))
        .thenThrow(new MakeFlightException("test exception"));

    assertThrows(
        InternalStairwayException.class, () -> jobService.enumerateJobs(userId, 10, null, null));
  }

  @Test
  void enumerateJobsInterruptedException() throws InterruptedException {
    String userId = "testUserId";

    when(mockStairway.getFlights((String) eq(null), eq(10), any(FlightFilter.class)))
        .thenThrow(new InterruptedException());

    // InterruptedException should be caught and re-thrown as an InternalStairwayException
    assertThrows(
        InternalStairwayException.class, () -> jobService.enumerateJobs(userId, 10, null, null));
  }

  @Test
  void retrieveAsyncJobResultRunning() throws InterruptedException {
    UUID jobId = TestUtils.TEST_NEW_UUID;
    FlightMap inputParameters = StairwayTestUtils.CREATE_JOB_INPUT_PARAMS;
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.RUNNING, jobId, inputParameters, new FlightMap());

    when(mockStairway.getFlightState(jobId.toString())).thenReturn(flightState);

    JobApiUtils.AsyncJobResult<String> result =
        jobService.retrieveAsyncJobResult(jobId, TestUtils.TEST_USER_ID_1, String.class, null);

    assertEquals(jobId.toString(), result.getJobReport().getId());
    assertEquals(202, result.getJobReport().getStatusCode());
    assertNull(result.getResult());
    assertNull(result.getApiErrorReport());
  }

  @Test
  void retrieveAsyncJobResultSucceeded() throws InterruptedException {
    UUID jobId = TestUtils.TEST_NEW_UUID;
    FlightMap inputParameters = StairwayTestUtils.CREATE_JOB_INPUT_PARAMS;
    FlightMap workingMap = new FlightMap();
    String testResponse = "test response";
    workingMap.put(JobMapKeys.RESPONSE, testResponse);
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.SUCCESS, jobId, inputParameters, workingMap);

    when(mockStairway.getFlightState(jobId.toString())).thenReturn(flightState);

    JobApiUtils.AsyncJobResult<String> result =
        jobService.retrieveAsyncJobResult(jobId, TestUtils.TEST_USER_ID_1, String.class, null);

    assertEquals(jobId.toString(), result.getJobReport().getId());
    assertEquals(200, result.getJobReport().getStatusCode());
    assertEquals(testResponse, result.getResult());
    assertNull(result.getApiErrorReport());
  }

  @Test
  void retrieveAsyncJobResultFailed() throws InterruptedException {
    UUID jobId = TestUtils.TEST_NEW_UUID;
    FlightMap inputParameters = StairwayTestUtils.CREATE_JOB_INPUT_PARAMS;
    // even on a fatal failure the response might have been written to the working map
    FlightMap workingMap = new FlightMap();
    String testResponse = "test response";
    workingMap.put(JobMapKeys.RESPONSE, testResponse);
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.ERROR, jobId, inputParameters, workingMap);
    String testErrorMsg = "test exception";
    RuntimeException exception = new RuntimeException(testErrorMsg);
    flightState.setException(exception);

    when(mockStairway.getFlightState(jobId.toString())).thenReturn(flightState);

    JobApiUtils.AsyncJobResult<String> result =
        jobService.retrieveAsyncJobResult(jobId, TestUtils.TEST_USER_ID_1, String.class, null);

    assertEquals(jobId.toString(), result.getJobReport().getId());
    assertEquals(500, result.getJobReport().getStatusCode());
    assertNull(result.getResult());
    assertEquals(testErrorMsg, result.getApiErrorReport().getMessage());
  }

  @Test
  void retrieveAsyncJobResultStairwayException() throws InterruptedException {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    // a MakeFlightException is an instance of StairwayException
    when(mockStairway.getFlightState(flightId.toString()))
        .thenThrow(new MakeFlightException("test exception"));

    assertThrows(
        InternalStairwayException.class,
        () ->
            jobService.retrieveAsyncJobResult(
                flightId, TestUtils.TEST_USER_ID_1, String.class, null));
  }
}
