package bio.terra.pipelines.controller;

import static bio.terra.pipelines.app.controller.JobApiUtils.*;
import static org.junit.jupiter.api.Assertions.*;

import bio.terra.common.exception.ErrorReportException;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.exception.InternalStairwayException;
import bio.terra.pipelines.dependencies.stairway.exception.InvalidResultStateException;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.generated.model.ApiErrorReport;
import bio.terra.pipelines.generated.model.ApiGetJobsResponse;
import bio.terra.pipelines.generated.model.ApiJobReport;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.springframework.http.HttpStatus;

class JobApiUtilsTest {

  @Test
  void testMapEnumeratedJobsToApi() {
    EnumeratedJobs bothJobs = StairwayTestUtils.ENUMERATED_JOBS;

    ApiGetJobsResponse mappedResponse = mapEnumeratedJobsToApi(bothJobs);

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

    ApiJobReport apiJobReport = mapFlightStateToApiJobReport(flightState);

    assertEquals(TestUtils.TEST_NEW_UUID.toString(), apiJobReport.getId());
    assertEquals(StairwayTestUtils.TEST_DESCRIPTION, apiJobReport.getDescription());
    assertEquals(
        "SUCCEEDED", apiJobReport.getStatus().name()); // "SUCCESS" gets mapped to "SUCCEEDED"
    assertEquals(StairwayTestUtils.TIME_SUBMITTED_1.toString(), apiJobReport.getSubmitted());
    assertEquals(StairwayTestUtils.TIME_COMPLETED_1.toString(), apiJobReport.getCompleted());
    assertEquals(TestUtils.TEST_RESULT_URL, apiJobReport.getResultURL());
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
        InvalidResultStateException.class, () -> mapFlightStateToApiJobReport(flightState));
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

    ApiJobReport apiJobReport = mapFlightStateToApiJobReport(flightState);

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

    ApiJobReport apiJobReport = mapFlightStateToApiJobReport(flightState);

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
        InvalidResultStateException.class, () -> mapFlightStateToApiJobReport(flightState));
  }

  @Test
  void testMapFlightStateToApiJobReportRunning() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.RUNNING, TestUtils.TEST_NEW_UUID);

    ApiJobReport apiJobReport = mapFlightStateToApiJobReport(flightState);

    assertEquals(TestUtils.TEST_NEW_UUID.toString(), apiJobReport.getId());
    assertEquals(StairwayTestUtils.TEST_DESCRIPTION, apiJobReport.getDescription());
    assertEquals("RUNNING", apiJobReport.getStatus().name());
    assertNull(apiJobReport.getCompleted());
    assertEquals(TestUtils.TEST_RESULT_URL, apiJobReport.getResultURL());
    assertEquals(202, apiJobReport.getStatusCode());
  }

  // the following tests are effectively tests of mapFlightStatusToApi(), which is private

  @Test
  void testMapFlightStateToApiJobReportRunningQueued() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.QUEUED, TestUtils.TEST_NEW_UUID);

    ApiJobReport apiJobReport = mapFlightStateToApiJobReport(flightState);

    assertEquals("RUNNING", apiJobReport.getStatus().name());
    assertNull(apiJobReport.getCompleted());
    assertEquals(202, apiJobReport.getStatusCode());
  }

  @Test
  void testMapFlightStateToApiJobReportRunningWaiting() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.WAITING, TestUtils.TEST_NEW_UUID);

    ApiJobReport apiJobReport = mapFlightStateToApiJobReport(flightState);

    assertEquals("RUNNING", apiJobReport.getStatus().name());
    assertNull(apiJobReport.getCompleted());
    assertEquals(202, apiJobReport.getStatusCode());
  }

  @Test
  void testMapFlightStateToApiJobReportRunningReady() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.READY, TestUtils.TEST_NEW_UUID);

    ApiJobReport apiJobReport = mapFlightStateToApiJobReport(flightState);

    assertEquals("RUNNING", apiJobReport.getStatus().name());
    assertNull(apiJobReport.getCompleted());
    assertEquals(202, apiJobReport.getStatusCode());
  }

  @Test
  void testMapFlightStateToApiJobReportRunningReadyToRestart() {
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.READY_TO_RESTART, TestUtils.TEST_NEW_UUID);

    ApiJobReport apiJobReport = mapFlightStateToApiJobReport(flightState);

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

    ApiJobReport apiJobReport = mapFlightStateToApiJobReport(flightState);

    assertEquals("FAILED", apiJobReport.getStatus().name());
    assertEquals(500, apiJobReport.getStatusCode());
  }

  @Test
  void testBuildApiErrorReportErrorReportExceptionCustomStatus() {
    String errorMessage = "some message";
    ErrorReportException exception =
        new ErrorReportException(errorMessage, null, HttpStatus.I_AM_A_TEAPOT) {};

    ApiErrorReport apiErrorReport = buildApiErrorReport(exception);

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

    ApiErrorReport apiErrorReport = buildApiErrorReport(exception);

    assertEquals(errorMessage, apiErrorReport.getMessage());
    assertEquals(418, apiErrorReport.getStatusCode());
    assertEquals(causes, apiErrorReport.getCauses());
  }

  @Test
  void testBuildApiErrorReport() {
    String errorMessage = "some message";
    InternalStairwayException exception = new InternalStairwayException(errorMessage);

    ApiErrorReport apiErrorReport = buildApiErrorReport(exception);

    assertEquals(errorMessage, apiErrorReport.getMessage());
    assertEquals(500, apiErrorReport.getStatusCode());
    assertTrue(apiErrorReport.getCauses().isEmpty());
  }
}
