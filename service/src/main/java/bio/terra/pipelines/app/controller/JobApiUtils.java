package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.ErrorReportException;
import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.dependencies.stairway.exception.InternalStairwayException;
import bio.terra.pipelines.dependencies.stairway.exception.InvalidResultStateException;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJob;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.generated.model.ApiErrorReport;
import bio.terra.pipelines.generated.model.ApiGetJobsResponse;
import bio.terra.pipelines.generated.model.ApiJobReport;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import bio.terra.stairway.exception.StairwayException;
import com.fasterxml.jackson.core.type.TypeReference;
import java.nio.file.Path;
import java.time.Instant;
import java.util.ArrayList;
import java.util.List;
import java.util.UUID;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.stereotype.Component;

@Component
public class JobApiUtils {
  private final JobService jobService;
  private final IngressConfiguration ingressConfiguration;

  @Autowired
  public JobApiUtils(JobService jobService, IngressConfiguration ingressConfiguration) {
    this.jobService = jobService;
    this.ingressConfiguration = ingressConfiguration;
  }

  public ApiGetJobsResponse mapEnumeratedJobsToApi(EnumeratedJobs enumeratedJobs) {
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

  public ApiJobReport mapFlightStateToApiJobReport(FlightState flightState) {
    FlightMap inputParameters = flightState.getInputParameters();
    String description = inputParameters.get(JobMapKeys.DESCRIPTION.getKeyName(), String.class);
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
          statusCode = resultMap.get(JobMapKeys.STATUS_CODE.getKeyName(), HttpStatus.class);
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

  public ApiErrorReport buildApiErrorReport(Exception exception) {
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
   * Retrieves the result of an asynchronous job.
   *
   * <p>Stairway has no concept of synchronous vs asynchronous flights. However, MC Terra has a
   * service-level standard result for asynchronous jobs which includes a ApiJobReport and either a
   * result or error if the job is complete. This is a convenience for callers who would otherwise
   * need to construct their own AsyncJobResult object.
   *
   * <p>Unlike retrieveJobResult, this will not throw for a flight in progress. Instead, it will
   * return a ApiJobReport without a result or error.
   */
  public <T> AsyncJobResult<T> retrieveAsyncJobResult(
      UUID jobId,
      String userId,
      PipelinesEnum pipelineId,
      Class<T> resultClass,
      TypeReference<T> typeReference) {
    try {
      FlightState flightState = jobService.retrieveJob(jobId, userId, pipelineId);
      ApiJobReport jobReport = mapFlightStateToApiJobReport(flightState);
      if (jobReport.getStatus().equals(ApiJobReport.StatusEnum.RUNNING)) {
        return new AsyncJobResult<T>().jobReport(jobReport);
      }

      // Job is complete, get the result
      JobService.JobResultOrException<T> resultOrException =
          jobService.retrieveJobResult(jobId, resultClass, typeReference);
      final ApiErrorReport errorReport;
      if (jobReport.getStatus().equals(ApiJobReport.StatusEnum.FAILED)) {
        errorReport = buildApiErrorReport(resultOrException.getException());
      } else {
        errorReport = null;
      }
      return new AsyncJobResult<T>()
          .jobReport(jobReport)
          .result(resultOrException.getResult())
          .errorReport(errorReport);
    } catch (StairwayException stairwayEx) {
      throw new InternalStairwayException(stairwayEx);
    }
  }

  private String resultUrlFromFlightState(FlightState flightState) {
    String resultPath =
        flightState.getInputParameters().get(JobMapKeys.RESULT_PATH.getKeyName(), String.class);
    if (resultPath == null) {
      resultPath = "";
    }
    // This is a little hacky, but GCP rejects non-https traffic and a local server does not
    // support it.
    String domainName = ingressConfiguration.getDomainName();
    String protocol = domainName.startsWith("localhost") ? "http://" : "https://";
    return protocol + Path.of(domainName, resultPath);
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
