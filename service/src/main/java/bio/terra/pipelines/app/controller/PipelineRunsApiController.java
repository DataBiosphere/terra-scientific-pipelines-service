package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.NotFoundException;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
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
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
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

  @Autowired
  public PipelineRunsApiController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      JobService jobService,
      PipelinesService pipelinesService,
      PipelineRunsService pipelineRunsService,
      PipelineInputsOutputsService pipelineInputsOutputsService,
      IngressConfiguration ingressConfiguration) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.jobService = jobService;
    this.pipelinesService = pipelinesService;
    this.pipelineRunsService = pipelineRunsService;
    this.pipelineInputsOutputsService = pipelineInputsOutputsService;
    this.ingressConfiguration = ingressConfiguration;
  }

  private static final Logger logger = LoggerFactory.getLogger(PipelineRunsApiController.class);

  private SamUser getAuthenticatedInfo() {
    return samUserFactory.from(request, samConfiguration.baseUri());
  }

  // PipelineRuns

  @Override
  public ResponseEntity<ApiPreparePipelineRunResponse> preparePipelineRun(
      @PathVariable("pipelineName") String pipelineName,
      @RequestBody ApiPreparePipelineRunRequestBody body) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    UUID jobId = body.getJobId();

    int pipelineVersion = body.getPipelineVersion();
    Map<String, Object> userProvidedInputs = body.getPipelineInputs();

    // validate the pipeline name and user-provided inputs
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);
    Pipeline pipeline = pipelinesService.getPipeline(validatedPipelineName);

    pipelinesService.validateUserProvidedInputs(
        pipeline.getPipelineInputDefinitions(), userProvidedInputs);

    logger.info(
        "Preparing {} pipeline (version {}) job (id {}) for user {} with validated inputs {}",
        pipelineName,
        pipelineVersion,
        jobId,
        userId,
        userProvidedInputs);

    Map<String, Map<String, String>> fileInputUploadUrls =
        pipelineRunsService.preparePipelineRun(pipeline, jobId, userId, userProvidedInputs);

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
   * @param pipelineName the pipeline to run
   * @param body the inputs for the pipeline
   * @return the created job response, which includes a job report containing the job ID,
   *     description, status, status code, submitted timestamp, completed timestamp (if completed),
   *     and result URL. The response also includes an error report if the job failed.
   */
  @Override
  public ResponseEntity<ApiAsyncPipelineRunResponse> startPipelineRun(
      @PathVariable("pipelineName") String pipelineName,
      @RequestBody ApiStartPipelineRunRequestBody body) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    UUID jobId = body.getJobControl().getId();

    String description = body.getDescription();

    // validate the pipeline name and user-provided inputs
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);
    Pipeline pipeline = pipelinesService.getPipeline(validatedPipelineName);

    logger.info("Starting {} pipeline job (id {}) for user {}", pipelineName, jobId, userId);

    PipelineRun pipelineRun =
        pipelineRunsService.startPipelineRun(pipeline, jobId, userId, description);

    ApiAsyncPipelineRunResponse createdRunResponse = pipelineRunToApi(pipelineRun, pipeline);

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
          "Pipeline run %s is still preparing, it has to be started before you can query the result"
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
  @Override
  public ResponseEntity<ApiGetPipelineRunsResponse> getAllPipelineRuns(
      Integer limit, String pageToken) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    int maxLimit = Math.min(limit, 100);

    // grab results from current page based on user provided inputs
    PageResponse<List<PipelineRun>> pageResults =
        pipelineRunsService.findPipelineRunsPaginated(maxLimit, pageToken, userId);

    // convert PageResponse object to list of ApiPipelineRun objects for response
    List<ApiPipelineRun> apiPipelineRuns =
        pageResults.content().stream()
            .map(
                pipelineRun ->
                    new ApiPipelineRun()
                        .jobId(pipelineRun.getJobId())
                        .status(pipelineRun.getStatus().name())
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
            .pageToken(pageResults.nextPageCursor());
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
            .wdlMethodVersion(
                pipelineRun
                    .getWdlMethodVersion())); // wdlMethodVersion comes from pipelineRun, since the
    // pipeline might have been updated since the pipelineRun began

    // if the pipeline run is successful, return the job report and add outputs to the response
    if (pipelineRun.getStatus().isSuccess()) {
      return response
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
                  .outputs(pipelineInputsOutputsService.formatPipelineRunOutputs(pipelineRun)));
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
