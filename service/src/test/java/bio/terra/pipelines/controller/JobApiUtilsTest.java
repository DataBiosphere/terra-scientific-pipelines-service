package bio.terra.pipelines.controller;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.ErrorReportException;
import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.app.controller.JobApiUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.dependencies.stairway.exception.InternalStairwayException;
import bio.terra.pipelines.dependencies.stairway.exception.InvalidResultStateException;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.generated.model.ApiErrorReport;
import bio.terra.pipelines.generated.model.ApiGetJobsResponse;
import bio.terra.pipelines.generated.model.ApiJobReport;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import bio.terra.stairway.exception.MakeFlightException;
import java.util.List;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.springframework.http.HttpStatus;

class JobApiUtilsTest extends BaseEmbeddedDbTest {

  @InjectMocks private JobApiUtils jobApiUtils;
  @Mock private JobService jobService;
  @Mock private IngressConfiguration ingressConfiguration;

  // we mock IngressConfiguration.getDomainName() to return "localhost", which maps to "http://"
  private final String fullResultURL =
      String.format("http://localhost/%s", TestUtils.TEST_RESULT_PATH);

  @BeforeEach
  void initMocks() {
    when(ingressConfiguration.getDomainName()).thenReturn("localhost");
  }

  @Test
  void testMapEnumeratedJobsToApi() {
    EnumeratedJobs bothJobs = StairwayTestUtils.ENUMERATED_JOBS;

    ApiGetJobsResponse mappedResponse = jobApiUtils.mapEnumeratedJobsToApi(bothJobs);

    assertEquals(bothJobs.getTotalResults(), mappedResponse.getTotalResults());
    assertEquals(bothJobs.getPageToken(), mappedResponse.getPageToken());
    assertEquals(bothJobs.getResults().size(), mappedResponse.getResults().size());
    assertEquals(
        bothJobs.getResults().get(0).getFlightState().getFlightId(),
        mappedResponse.getResults().get(0).getId());
    assertEquals(
        bothJobs.getResults().get(1).getFlightState().getFlightId(),
        mappedResponse.getResults().get(1).getId());
  }

  @Test
  void testMapFlightStateToApiJobReportSucceeded() {
    FlightState flightState = StairwayTestUtils.FLIGHT_STATE_DONE_SUCCESS_1;

    ApiJobReport apiJobReport = jobApiUtils.mapFlightStateToApiJobReport(flightState);

    assertEquals(TestUtils.TEST_NEW_UUID.toString(), apiJobReport.getId());
    assertEquals(StairwayTestUtils.TEST_DESCRIPTION, apiJobReport.getDescription());
    assertEquals(
        "SUCCEEDED", apiJobReport.getStatus().name()); // "SUCCESS" gets mapped to "SUCCEEDED"
    assertEquals(StairwayTestUtils.TIME_SUBMITTED_1.toString(), apiJobReport.getSubmitted());
    assertEquals(StairwayTestUtils.TIME_COMPLETED_1.toString(), apiJobReport.getCompleted());
    assertEquals(fullResultURL, apiJobReport.getResultURL());
    // if there is no status code in the working map, we assume it's a success/200
    assertEquals(200, apiJobReport.getStatusCode());
  }

  @Test
  void testMapFlightStateToApiJobReportSucceededNoCompletedTime() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.SUCCESS,
            TestUtils.TEST_NEW_UUID,
            StairwayTestUtils.CREATE_JOB_INPUT_PARAMS,
            StairwayTestUtils.EMPTY_WORKING_MAP,
            StairwayTestUtils.TIME_SUBMITTED_1,
            null);

    // a completed job should have a completed time, otherwise it's an error
    assertThrows(
        InvalidResultStateException.class,
        () -> jobApiUtils.mapFlightStateToApiJobReport(flightState));
  }

  @Test
  void testMapFlightStateToApiJobReportSucceededWithStatusCode() {
    // Ensure the custom status code in the working map gets extracted and used
    HttpStatus httpStatus = HttpStatus.I_AM_A_TEAPOT; // status code 418
    FlightMap flightMapWithStatusCode = new FlightMap();
    flightMapWithStatusCode.put(JobMapKeys.STATUS_CODE.getKeyName(), httpStatus);

    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.SUCCESS,
            TestUtils.TEST_NEW_UUID,
            StairwayTestUtils.CREATE_JOB_INPUT_PARAMS,
            flightMapWithStatusCode,
            StairwayTestUtils.TIME_SUBMITTED_1,
            StairwayTestUtils.TIME_COMPLETED_1);

    ApiJobReport apiJobReport = jobApiUtils.mapFlightStateToApiJobReport(flightState);

    assertEquals(
        "SUCCEEDED", apiJobReport.getStatus().name()); // "SUCCESS" gets mapped to "SUCCEEDED"
    assertEquals(418, apiJobReport.getStatusCode());
  }

  @Test
  void testMapFlightStateToApiJobReportFailed() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.ERROR, TestUtils.TEST_NEW_UUID);
    flightState.setException(new InternalStairwayException("some message"));

    ApiJobReport apiJobReport = jobApiUtils.mapFlightStateToApiJobReport(flightState);

    assertEquals(TestUtils.TEST_NEW_UUID.toString(), apiJobReport.getId());
    assertEquals("FAILED", apiJobReport.getStatus().name()); // ERROR gets mapped to FAILED
    // InternalServerError is a 500
    assertEquals(500, apiJobReport.getStatusCode());
  }

  @Test
  void testMapFlightStateToApiJobReportFailedMissingException() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.ERROR,
            TestUtils.TEST_NEW_UUID,
            StairwayTestUtils.CREATE_JOB_INPUT_PARAMS,
            StairwayTestUtils.EMPTY_WORKING_MAP,
            StairwayTestUtils.TIME_SUBMITTED_1,
            StairwayTestUtils.TIME_COMPLETED_1);

    assertThrows(
        InvalidResultStateException.class,
        () -> jobApiUtils.mapFlightStateToApiJobReport(flightState));
  }

  @Test
  void testMapFlightStateToApiJobReportRunning() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.RUNNING, TestUtils.TEST_NEW_UUID);

    ApiJobReport apiJobReport = jobApiUtils.mapFlightStateToApiJobReport(flightState);

    assertEquals(TestUtils.TEST_NEW_UUID.toString(), apiJobReport.getId());
    assertEquals(StairwayTestUtils.TEST_DESCRIPTION, apiJobReport.getDescription());
    assertEquals("RUNNING", apiJobReport.getStatus().name());
    assertNull(apiJobReport.getCompleted());
    assertEquals(fullResultURL, apiJobReport.getResultURL());
    assertEquals(202, apiJobReport.getStatusCode());
  }

  // the following tests are effectively tests of mapFlightStatusToApi(), which is private

  @Test
  void testMapFlightStateToApiJobReportRunningQueued() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.QUEUED, TestUtils.TEST_NEW_UUID);

    ApiJobReport apiJobReport = jobApiUtils.mapFlightStateToApiJobReport(flightState);

    assertEquals("RUNNING", apiJobReport.getStatus().name());
    assertNull(apiJobReport.getCompleted());
    assertEquals(202, apiJobReport.getStatusCode());
  }

  @Test
  void testMapFlightStateToApiJobReportRunningWaiting() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.WAITING, TestUtils.TEST_NEW_UUID);

    ApiJobReport apiJobReport = jobApiUtils.mapFlightStateToApiJobReport(flightState);

    assertEquals("RUNNING", apiJobReport.getStatus().name());
    assertNull(apiJobReport.getCompleted());
    assertEquals(202, apiJobReport.getStatusCode());
  }

  @Test
  void testMapFlightStateToApiJobReportRunningReady() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.READY, TestUtils.TEST_NEW_UUID);

    ApiJobReport apiJobReport = jobApiUtils.mapFlightStateToApiJobReport(flightState);

    assertEquals("RUNNING", apiJobReport.getStatus().name());
    assertNull(apiJobReport.getCompleted());
    assertEquals(202, apiJobReport.getStatusCode());
  }

  @Test
  void testMapFlightStateToApiJobReportRunningReadyToRestart() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.READY_TO_RESTART, TestUtils.TEST_NEW_UUID);

    ApiJobReport apiJobReport = jobApiUtils.mapFlightStateToApiJobReport(flightState);

    assertEquals("RUNNING", apiJobReport.getStatus().name());
    assertNull(apiJobReport.getCompleted());
    assertEquals(202, apiJobReport.getStatusCode());
  }

  @Test
  void testMapFlightStateToApiJobReportFailedFatal() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.FATAL, TestUtils.TEST_NEW_UUID);
    flightState.setException(new InternalStairwayException("some message"));

    ApiJobReport apiJobReport = jobApiUtils.mapFlightStateToApiJobReport(flightState);

    assertEquals("FAILED", apiJobReport.getStatus().name());
    assertEquals(500, apiJobReport.getStatusCode());
  }

  @Test
  void testBuildApiErrorReportErrorReportExceptionCustomStatus() {
    String errorMessage = "some message";
    ErrorReportException exception =
        new ErrorReportException(errorMessage, null, HttpStatus.I_AM_A_TEAPOT) {};

    ApiErrorReport apiErrorReport = jobApiUtils.buildApiErrorReport(exception);

    assertEquals(errorMessage, apiErrorReport.getMessage());
    assertEquals(418, apiErrorReport.getStatusCode());
    assertTrue(apiErrorReport.getCauses().isEmpty());
  }

  @Test
  void testBuildApiErrorReportErrorReportExceptionCustomStatusWithCauses() {
    String errorMessage = "some message";
    List<String> causes = List.of("cause 1", "cause 2");
    ErrorReportException exception =
        new ErrorReportException(errorMessage, causes, HttpStatus.I_AM_A_TEAPOT) {};

    ApiErrorReport apiErrorReport = jobApiUtils.buildApiErrorReport(exception);

    assertEquals(errorMessage, apiErrorReport.getMessage());
    assertEquals(418, apiErrorReport.getStatusCode());
    assertEquals(causes, apiErrorReport.getCauses());
  }

  @Test
  void testBuildApiErrorReport() {
    String errorMessage = "some message";
    InternalStairwayException exception = new InternalStairwayException(errorMessage);

    ApiErrorReport apiErrorReport = jobApiUtils.buildApiErrorReport(exception);

    assertEquals(errorMessage, apiErrorReport.getMessage());
    assertEquals(500, apiErrorReport.getStatusCode());
    assertTrue(apiErrorReport.getCauses().isEmpty());
  }

  @Test
  void testRetrieveAsyncJobResultRunning() throws InterruptedException {
    UUID jobId = TestUtils.TEST_NEW_UUID;
    FlightMap inputParameters = StairwayTestUtils.CREATE_JOB_INPUT_PARAMS;
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.RUNNING, jobId, inputParameters, new FlightMap());

    when(jobService.retrieveJob(any(), any(), any())).thenReturn(flightState);

    JobApiUtils.AsyncJobResult<String> result =
        jobApiUtils.retrieveAsyncJobResult(
            jobId, TestUtils.TEST_USER_ID_1, PipelinesEnum.IMPUTATION_MINIMAC4, String.class, null);

    assertEquals(jobId.toString(), result.getJobReport().getId());
    assertEquals(202, result.getJobReport().getStatusCode());
    assertNull(result.getResult());
    assertNull(result.getApiErrorReport());
  }

  @Test
  void testRetrieveAsyncJobResultSucceeded() {
    UUID jobId = TestUtils.TEST_NEW_UUID;
    FlightMap inputParameters = StairwayTestUtils.CREATE_JOB_INPUT_PARAMS;
    FlightMap workingMap = new FlightMap();
    String testResponse = "test response";
    workingMap.put(JobMapKeys.RESPONSE.getKeyName(), testResponse);
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.SUCCESS, jobId, inputParameters, workingMap);

    when(jobService.retrieveJob(any(), any(), any())).thenReturn(flightState);
    when(jobService.retrieveJobResult(any(), any(), any()))
        .thenReturn(new JobService.JobResultOrException<>().result(testResponse));

    JobApiUtils.AsyncJobResult<String> result =
        jobApiUtils.retrieveAsyncJobResult(
            jobId, TestUtils.TEST_USER_ID_1, PipelinesEnum.IMPUTATION_MINIMAC4, String.class, null);

    assertEquals(jobId.toString(), result.getJobReport().getId());
    assertEquals(200, result.getJobReport().getStatusCode());
    assertEquals(testResponse, result.getResult());
    assertNull(result.getApiErrorReport());
  }

  @Test
  void testRetrieveAsyncJobResultFailed() {
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
    RuntimeException exception = new RuntimeException(testErrorMsg);
    flightState.setException(exception);

    when(jobService.retrieveJob(any(), any(), any())).thenReturn(flightState);
    when(jobService.retrieveJobResult(any(), any(), any()))
        .thenReturn(new JobService.JobResultOrException<>().exception(exception));

    JobApiUtils.AsyncJobResult<String> result =
        jobApiUtils.retrieveAsyncJobResult(
            jobId, TestUtils.TEST_USER_ID_1, PipelinesEnum.IMPUTATION_MINIMAC4, String.class, null);

    assertEquals(jobId.toString(), result.getJobReport().getId());
    assertEquals(500, result.getJobReport().getStatusCode());
    assertNull(result.getResult());
    assertEquals(testErrorMsg, result.getApiErrorReport().getMessage());
  }

  @Test
  void retrieveAsyncJobResultStairwayException() {
    UUID flightId = TestUtils.TEST_NEW_UUID;
    // a MakeFlightException is an instance of StairwayException
    when(jobService.retrieveJob(any(), any(), any()))
        .thenThrow(new MakeFlightException("test exception"));

    assertThrows(
        InternalStairwayException.class,
        () ->
            jobApiUtils.retrieveAsyncJobResult(
                flightId,
                TestUtils.TEST_USER_ID_1,
                TestUtils.TEST_PIPELINE_1_ENUM,
                String.class,
                null));
  }
}
