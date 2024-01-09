package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.ErrorReportException;
import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.dependencies.stairway.StairwayJobMapKeys;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.dependencies.stairway.exception.InvalidResultStateException;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJob;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.generated.model.ApiErrorReport;
import bio.terra.pipelines.generated.model.ApiGetJobsResponse;
import bio.terra.pipelines.generated.model.ApiJobReport;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import java.time.Instant;
import java.util.ArrayList;
import java.util.List;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.stereotype.Component;

@Component
public class JobApiUtils {
  private final StairwayJobService stairwayJobService;
  private final IngressConfiguration ingressConfig;

  @Autowired
  JobApiUtils(StairwayJobService stairwayJobService, IngressConfiguration ingressConfig) {
    this.stairwayJobService = stairwayJobService;
    this.ingressConfig = ingressConfig;
  }

  public static ApiGetJobsResponse mapEnumeratedJobsToApi(EnumeratedJobs enumeratedJobs) {

    // Convert the result to API-speak
    List<ApiJobReport> apiJobList = new ArrayList<>();
    for (EnumeratedJob enumeratedJob : enumeratedJobs.getResults()) {
      ApiJobReport jobReport = mapFlightStateToApiJobReport(enumeratedJob.getFlightState());
      apiJobList.add(jobReport);
    }

    return new ApiGetJobsResponse()
        .pageToken(enumeratedJobs.getPageToken())
        .totalResults(enumeratedJobs.getTotalResults())
        .results(apiJobList);
  }

  public static ApiJobReport mapFlightStateToApiJobReport(FlightState flightState) {
    FlightMap inputParameters = flightState.getInputParameters();
    String description =
        inputParameters.get(StairwayJobMapKeys.DESCRIPTION.getKeyName(), String.class);
    FlightStatus flightStatus = flightState.getFlightStatus();
    String submittedDate = flightState.getSubmitted().toString();
    ApiJobReport.StatusEnum jobStatus = mapFlightStatusToApi(flightStatus);

    String completedDate = null;
    HttpStatus statusCode = HttpStatus.ACCEPTED;

    if (jobStatus != ApiJobReport.StatusEnum.RUNNING) {
      // If the job is completed, the JobReport should include a result code indicating success or
      // failure. For failed jobs, this code is the error code. For successful jobs, this is the
      // code specified by the flight if present, or a default of 200 if not.
      completedDate =
          flightState
              .getCompleted()
              .map(Instant::toString)
              .orElseThrow(
                  () -> new InvalidResultStateException("No completed time for completed flight"));
      switch (jobStatus) {
        case FAILED -> {
          int errorCode =
              flightState
                  .getException()
                  .map(e -> buildApiErrorReport(e).getStatusCode())
                  .orElseThrow(
                      () ->
                          new InvalidResultStateException(
                              String.format(
                                  "Flight %s failed with no exception reported",
                                  flightState.getFlightId())));
          statusCode = HttpStatus.valueOf(errorCode);
        }
        case SUCCEEDED -> {
          FlightMap resultMap =
              flightState.getResultMap().orElseThrow(InvalidResultStateException::noResultMap);
          statusCode = resultMap.get(StairwayJobMapKeys.STATUS_CODE.getKeyName(), HttpStatus.class);
          if (statusCode == null) {
            statusCode = HttpStatus.OK;
          }
        }
        default -> throw new IllegalStateException(
            "Cannot get status code of flight in unknown state " + jobStatus);
      }
    }

    return new ApiJobReport()
        .id(flightState.getFlightId())
        .description(description)
        .status(jobStatus)
        .statusCode(statusCode.value())
        .submitted(submittedDate)
        .completed(completedDate)
        .resultURL(resultUrlFromFlightState(flightState));
  }

  private static ApiJobReport.StatusEnum mapFlightStatusToApi(FlightStatus flightStatus) {
    switch (flightStatus) {
      case RUNNING:
      case QUEUED:
      case WAITING:
      case READY:
      case READY_TO_RESTART:
        return ApiJobReport.StatusEnum.RUNNING;
      case SUCCESS:
        return ApiJobReport.StatusEnum.SUCCEEDED;
      case ERROR:
      case FATAL:
      default:
        return ApiJobReport.StatusEnum.FAILED;
    }
  }

  public static ApiErrorReport buildApiErrorReport(Exception exception) {
    if (exception instanceof ErrorReportException errorReport) {
      return new ApiErrorReport()
          .message(errorReport.getMessage())
          .statusCode(errorReport.getStatusCode().value())
          .causes(errorReport.getCauses());
    } else {
      return new ApiErrorReport()
          .message(exception.getMessage())
          .statusCode(HttpStatus.INTERNAL_SERVER_ERROR.value())
          .causes(null);
    }
  }

  private static String resultUrlFromFlightState(FlightState flightState) {
    String resultPath =
        flightState
            .getInputParameters()
            .get(StairwayJobMapKeys.RESULT_PATH.getKeyName(), String.class);
    if (resultPath == null) {
      resultPath = "";
    }
    // This is a little hacky, but GCP rejects non-https traffic and a local server does not
    // support it.
    //    String protocol =
    //        ingressConfig.getDomainName().startsWith("localhost") ? "http://" : "https://";
    return resultPath;
    //    return protocol + Path.of(ingressConfig.getDomainName(), resultPath);
  }
}
