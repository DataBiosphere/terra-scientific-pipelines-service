package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.ErrorReportException;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.exception.InvalidResultStateException;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJob;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.generated.model.ApiErrorReport;
import bio.terra.pipelines.generated.model.ApiGetJobsResponse;
import bio.terra.pipelines.generated.model.ApiJobReport;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import java.nio.file.Path;
import java.time.Instant;
import java.util.ArrayList;
import java.util.List;
import java.util.UUID;
import org.springframework.http.HttpStatus;

public class JobApiUtils {
  private JobApiUtils() {}

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

  @SuppressWarnings(
      "java:S2259") // suppress warning for possible NPE when calling pipelineName.getValue(),
  //  since we do validate that pipelineName is not null in `validateRequiredEntries`
  public static ApiJobReport mapFlightStateToApiJobReport(FlightState flightState) {
    FlightMap inputParameters = flightState.getInputParameters();

    // Validate that the required entries are present in the input parameters
    FlightUtils.validateRequiredEntries(inputParameters, JobMapKeys.DOMAIN_NAME);
    String domainName = inputParameters.get(JobMapKeys.DOMAIN_NAME, String.class);

    // description is optional
    String description = inputParameters.get(JobMapKeys.DESCRIPTION, String.class);

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
          statusCode = resultMap.get(JobMapKeys.STATUS_CODE, HttpStatus.class);
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
        .resultURL(getAsyncResultEndpoint(domainName, UUID.fromString(flightState.getFlightId())));
  }

  private static ApiJobReport.StatusEnum mapFlightStatusToApi(FlightStatus flightStatus) {
    switch (flightStatus) {
      case RUNNING, QUEUED, WAITING, READY, READY_TO_RESTART:
        return ApiJobReport.StatusEnum.RUNNING;
      case SUCCESS:
        return ApiJobReport.StatusEnum.SUCCEEDED;
      case ERROR, FATAL:
        return ApiJobReport.StatusEnum.FAILED;
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

  /**
   * Returns the result endpoint corresponding to an async request. The endpoint is used to build an
   * ApiJobReport. This method generates a result endpoin with the form:
   * {protocol}{domainName}/api/pipelineruns/v1/result/{jobId}. The domainName should come from the
   * IngressConfiguration.
   *
   * @param domainName the domain name from the ingress configuration
   * @param jobId identifier for the job
   * @return a string with the result endpoint URL
   */
  public static String getAsyncResultEndpoint(String domainName, UUID jobId) {
    String endpointPath = "/api/pipelineruns/v1/result/%s".formatted(jobId);

    // This is a little hacky, but GCP rejects non-https traffic and a local server does not
    // support it.
    String protocol = domainName.startsWith("localhost") ? "http://" : "https://";

    return protocol + Path.of(domainName, endpointPath);
  }

  /**
   * The API result of an asynchronous job is a ApiJobReport and exactly one of a job result of an
   * ApiErrorReport. If the job is incomplete, only jobReport will be present.
   *
   * @param <T> Class of the result object
   */
  public static class AsyncJobResult<T> {
    private ApiJobReport jobReport;
    private T result;
    private ApiErrorReport errorReport;

    public T getResult() {
      return result;
    }

    public AsyncJobResult<T> result(T result) {
      this.result = result;
      return this;
    }

    public ApiErrorReport getApiErrorReport() {
      return errorReport;
    }

    public AsyncJobResult<T> errorReport(ApiErrorReport errorReport) {
      this.errorReport = errorReport;
      return this;
    }

    public ApiJobReport getJobReport() {
      return jobReport;
    }

    public AsyncJobResult<T> jobReport(ApiJobReport jobReport) {
      this.jobReport = jobReport;
      return this;
    }
  }
}
