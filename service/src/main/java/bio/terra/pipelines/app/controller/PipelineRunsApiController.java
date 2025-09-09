package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.NotFoundException;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.app.configuration.internal.PipelinesCommonConfiguration;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.pagination.PageResponse;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.generated.api.PipelineRunsApi;
import bio.terra.pipelines.generated.model.*;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.PipelinesService;
import io.swagger.annotations.Api;
import jakarta.servlet.http.HttpServletRequest;
import java.time.Instant;
import java.time.temporal.ChronoUnit;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import java.util.stream.Collectors;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.data.domain.Page;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestBody;

/** PipelineRuns controller */
@Controller
@Api(tags = {"pipelineRuns"})
public class PipelineRunsApiController implements PipelineRunsApi {
  private final SamConfiguration samConfiguration;
  private final SamUserFactory samUserFactory;
  private final HttpServletRequest request;
  private final JobService jobService;
  private final PipelinesService pipelinesService;
  private final PipelineRunsService pipelineRunsService;
  private final PipelineInputsOutputsService pipelineInputsOutputsService;
  private final IngressConfiguration ingressConfiguration;
  private final PipelinesCommonConfiguration pipelinesCommonConfiguration;

  @Autowired
  public PipelineRunsApiController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      JobService jobService,
      PipelinesService pipelinesService,
      PipelineRunsService pipelineRunsService,
      PipelineInputsOutputsService pipelineInputsOutputsService,
      IngressConfiguration ingressConfiguration,
      PipelinesCommonConfiguration pipelinesCommonConfiguration) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.jobService = jobService;
    this.pipelinesService = pipelinesService;
    this.pipelineRunsService = pipelineRunsService;
    this.pipelineInputsOutputsService = pipelineInputsOutputsService;
    this.ingressConfiguration = ingressConfiguration;
    this.pipelinesCommonConfiguration = pipelinesCommonConfiguration;
  }

  private static final Logger logger = LoggerFactory.getLogger(PipelineRunsApiController.class);

  private SamUser getAuthenticatedInfo() {
    return samUserFactory.from(request, samConfiguration.baseUri());
  }

  // PipelineRuns

  @Override
  public ResponseEntity<ApiPreparePipelineRunResponse> preparePipelineRun(
      @RequestBody ApiPreparePipelineRunRequestBody body) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    UUID jobId = body.getJobId();
    String pipelineName = body.getPipelineName();
    String description = body.getDescription();

    Integer pipelineVersion = body.getPipelineVersion();
    Map<String, Object> userProvidedInputs = body.getPipelineInputs();

    // validate the pipeline name and user-provided inputs
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);
    Pipeline pipeline = pipelinesService.getPipeline(validatedPipelineName, pipelineVersion);

    pipelineInputsOutputsService.validateUserProvidedInputs(
        pipeline.getPipelineInputDefinitions(), userProvidedInputs);

    logger.info(
        "Preparing {} pipeline (version {}) job (id {}) for user {} with validated inputs {}",
        pipelineName,
        pipelineVersion,
        jobId,
        userId,
        userProvidedInputs);

    Map<String, Map<String, String>> fileInputUploadUrls =
        pipelineRunsService.preparePipelineRun(
            pipeline, jobId, userId, userProvidedInputs, description);

    ApiPreparePipelineRunResponse prepareResponse =
        new ApiPreparePipelineRunResponse().jobId(jobId).fileInputUploadUrls(fileInputUploadUrls);

    return new ResponseEntity<>(prepareResponse, HttpStatus.OK);
  }

  /**
   * Kicks off the asynchronous process (managed by Stairway) of gathering user-provided inputs,
   * running the specified pipeline, and delivering the outputs to the user.
   *
   * <p>The run is created with a user-provided job ID (uuid).
   *
   * @param body the inputs for the pipeline
   * @return the created job response, which includes a job report containing the job ID,
   *     description, status, status code, submitted timestamp, completed timestamp (if completed),
   *     and result URL. The response also includes an error report if the job failed.
   */
  @Override
  public ResponseEntity<ApiAsyncPipelineRunResponse> startPipelineRun(
      @RequestBody ApiStartPipelineRunRequestBody body) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    UUID jobId = body.getJobControl().getId();

    PipelineRun pipelineRunBeforeStart = pipelineRunsService.getPipelineRun(jobId, userId);
    Pipeline pipeline = pipelinesService.getPipelineById(pipelineRunBeforeStart.getPipelineId());

    logger.info(
        "Starting {} pipeline job (id {}) for user {}",
        pipeline.getName().getValue(),
        jobId,
        userId);

    // we don't want to start a pipeline whose input data would have been deleted because
    // of the TTL on user data.
    if (pipelineRunBeforeStart
        .getCreated()
        .plus(pipelinesCommonConfiguration.getUserDataTtlDays(), ChronoUnit.DAYS)
        .isBefore(Instant.now())) {
      throw new BadRequestException(
          "Pipeline run was prepared more than %s days ago; it cannot be started"
              .formatted(pipelinesCommonConfiguration.getUserDataTtlDays()));
    }

    PipelineRun pipelineRunAfterStart =
        pipelineRunsService.startPipelineRun(pipeline, jobId, userId);

    ApiAsyncPipelineRunResponse createdRunResponse =
        pipelineRunToApi(pipelineRunAfterStart, pipeline);

    return new ResponseEntity<>(
        createdRunResponse, getAsyncResponseCode(createdRunResponse.getJobReport()));
  }

  @Override
  public ResponseEntity<ApiAsyncPipelineRunResponse> getPipelineRunResult(
      @PathVariable("jobId") UUID jobId) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();

    PipelineRun pipelineRun = pipelineRunsService.getPipelineRun(jobId, userId);
    if (pipelineRun == null) {
      throw new NotFoundException("Pipeline run %s not found".formatted(jobId));
    }

    if (pipelineRun.getStatus().equals(CommonPipelineRunStatusEnum.PREPARING)) {
      throw new BadRequestException(
          "Pipeline run %s is still preparing; it has to be started before you can query the result"
              .formatted(jobId));
    }

    Pipeline pipeline = pipelinesService.getPipelineById(pipelineRun.getPipelineId());

    ApiAsyncPipelineRunResponse runResponse = pipelineRunToApi(pipelineRun, pipeline);

    return new ResponseEntity<>(runResponse, getAsyncResponseCode(runResponse.getJobReport()));
  }

  /**
   * Returns a paginated list of Pipeline Runs for the user calling this endpoint. This will return
   * in descending order i.e. will return the most recent runs first. pageToken is used to navigate
   * to the corresponding "page" of results using <a
   * href="https://bun.uptrace.dev/guide/cursor-pagination.html">cursor based pagination</a>
   *
   * @param limit - how many results the caller wants to be returned, maximum value is 100
   * @param pageToken - token used for cursor based pagination. if not supplied then the first page
   *     will be returned
   * @return ResponseEntity containing the current page of results and a page token for the next
   *     page if a next page exists
   */
  @Deprecated
  @Override
  public ResponseEntity<ApiGetPipelineRunsResponse> getAllPipelineRuns(
      Integer limit, String pageToken) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    int maxLimit = Math.min(limit, 100);

    // grab results from current page based on user provided inputs
    PageResponse<List<PipelineRun>> pageResults =
        pipelineRunsService.findPipelineRunsPaginated(maxLimit, pageToken, userId);

    // convert list of pipelines to map of id to pipeline
    Map<Long, Pipeline> pipelineIdToPipeline =
        pipelinesService.getPipelines().stream().collect(Collectors.toMap(Pipeline::getId, p -> p));

    int totalResults = Math.toIntExact(pipelineRunsService.getPipelineRunCount(userId));

    // convert PageResponse object to list of ApiPipelineRun objects for response
    List<ApiPipelineRun> apiPipelineRuns =
        pageResults.content().stream()
            .map(
                pipelineRun ->
                    new ApiPipelineRun()
                        .jobId(pipelineRun.getJobId())
                        .pipelineName(
                            pipelineIdToPipeline
                                .get(pipelineRun.getPipelineId())
                                .getName()
                                .getValue())
                        .pipelineVersion(
                            pipelineIdToPipeline.get(pipelineRun.getPipelineId()).getVersion())
                        .status(pipelineRun.getStatus().name())
                        .quotaConsumed(pipelineRun.getQuotaConsumed())
                        .description(pipelineRun.getDescription())
                        .timeSubmitted(pipelineRun.getCreated().toString())
                        .timeCompleted(
                            pipelineRun.getStatus().isCompleted()
                                ? pipelineRun.getUpdated().toString()
                                : null))
            .toList();

    ApiGetPipelineRunsResponse apiGetPipelineRunsResponse =
        new ApiGetPipelineRunsResponse()
            .results(apiPipelineRuns)
            .totalResults(totalResults)
            .pageToken(pageResults.nextPageCursor());
    return new ResponseEntity<>(apiGetPipelineRunsResponse, HttpStatus.OK);
  }

  // TODO: check indices
  @Override
  public ResponseEntity<ApiGetPipelineRunsResponseV2> getAllPipelineRunsV2(
      Integer pageNumber, Integer pageSize, String sortProperty, String sortDirection) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    int maxPageSize = Math.min(pageSize, 100);

    // grab results from current page based on user provided inputs
    // PageRequest is zero-indexed, but for the API we want to be one-indexed for user-friendliness
    Page<PipelineRun> pageResults =
        pipelineRunsService.findPipelineRunsPaginated(
            pageNumber - 1, maxPageSize, sortProperty, sortDirection, userId);

    // convert list of pipelines to map of id to pipeline
    Map<Long, Pipeline> pipelineIdToPipeline =
        pipelinesService.getPipelines().stream().collect(Collectors.toMap(Pipeline::getId, p -> p));

    int totalResults = Math.toIntExact(pipelineRunsService.getPipelineRunCount(userId));

    // convert Page object to list of ApiPipelineRun objects for response
    List<ApiPipelineRun> apiPipelineRuns =
        pageResults
            .get()
            .map(
                pipelineRun ->
                    new ApiPipelineRun()
                        .jobId(pipelineRun.getJobId())
                        .pipelineName(
                            pipelineIdToPipeline
                                .get(pipelineRun.getPipelineId())
                                .getName()
                                .getValue())
                        .pipelineVersion(
                            pipelineIdToPipeline.get(pipelineRun.getPipelineId()).getVersion())
                        .status(pipelineRun.getStatus().name())
                        .quotaConsumed(pipelineRun.getQuotaConsumed())
                        .description(pipelineRun.getDescription())
                        .timeSubmitted(pipelineRun.getCreated().toString())
                        .timeCompleted(
                            pipelineRun.getStatus().isCompleted()
                                ? pipelineRun.getUpdated().toString()
                                : null))
            .toList();

    ApiGetPipelineRunsResponseV2 apiGetPipelineRunsResponse =
        new ApiGetPipelineRunsResponseV2().results(apiPipelineRuns).totalResults(totalResults);
    return new ResponseEntity<>(apiGetPipelineRunsResponse, HttpStatus.OK);
  }

  // helper methods
  /**
   * Converts a PipelineRun to an ApiAsyncPipelineRunResponse. If the PipelineRun has completed
   * successfully (isSuccess is true), we know there are no errors to retrieve and all the
   * information needed to return to the user is available in the pipeline_runs table.
   *
   * <p>If the PipelineRun is not marked as a success, we retrieve the running and/or error
   * information from Stairway.
   *
   * @param pipelineRun the PipelineRun to convert
   * @return ApiAsyncPipelineRunResponse
   */
  private ApiAsyncPipelineRunResponse pipelineRunToApi(PipelineRun pipelineRun, Pipeline pipeline) {
    ApiAsyncPipelineRunResponse response = new ApiAsyncPipelineRunResponse();
    response.pipelineRunReport(
        new ApiPipelineRunReport()
            .pipelineName(pipeline.getName().getValue())
            .pipelineVersion(pipeline.getVersion())
            .toolVersion(
                pipelineRun.getToolVersion())); // toolVersion comes from pipelineRun, since the
    // pipeline might have been updated since the pipelineRun began

    // if the pipeline run is successful, return the job report and add outputs to the response
    if (pipelineRun.getStatus().isSuccess()) {
      // calculate the expiration date for the output files
      Instant outputExpirationDate =
          pipelineRun
              .getUpdated()
              .plus(pipelinesCommonConfiguration.getUserDataTtlDays(), ChronoUnit.DAYS);
      response
          .jobReport(
              new ApiJobReport()
                  .id(pipelineRun.getJobId().toString())
                  .description(pipelineRun.getDescription())
                  .status(ApiJobReport.StatusEnum.SUCCEEDED)
                  .statusCode(HttpStatus.OK.value())
                  .submitted(pipelineRun.getCreated().toString())
                  .completed(pipelineRun.getUpdated().toString())
                  .resultURL(
                      JobApiUtils.getAsyncResultEndpoint(
                          ingressConfiguration.getDomainName(), pipelineRun.getJobId())))
          .pipelineRunReport(
              response
                  .getPipelineRunReport()
                  .outputExpirationDate(outputExpirationDate.toString()));

      // return outputs if we have not passed the output expiration date
      if (outputExpirationDate.isAfter(Instant.now())) {
        response
            .getPipelineRunReport()
            .outputs(pipelineInputsOutputsService.formatPipelineRunOutputs(pipelineRun));
      }
      return response;

    } else {
      JobApiUtils.AsyncJobResult<String> jobResult =
          jobService.retrieveAsyncJobResult(
              pipelineRun.getJobId(), pipelineRun.getUserId(), String.class, null);
      return response
          .jobReport(jobResult.getJobReport())
          .errorReport(jobResult.getApiErrorReport());
    }
  }

  /**
   * Return the appropriate response code for an endpoint, given an async job report. For a job
   * that's still running, this is 202. For a job that's finished (either succeeded or failed), the
   * endpoint should return 200. More informational status codes will be included in either the
   * response or error report bodies.
   */
  public static HttpStatus getAsyncResponseCode(ApiJobReport jobReport) {
    return jobReport.getStatus() == ApiJobReport.StatusEnum.RUNNING
        ? HttpStatus.ACCEPTED
        : HttpStatus.OK;
  }
}
