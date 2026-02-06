package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.NotFoundException;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelineRunFilterSpecification;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.generated.api.PipelineRunsApi;
import bio.terra.pipelines.generated.model.*;
import bio.terra.pipelines.service.DownloadCallCounterService;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.service.QuotasService;
import io.swagger.annotations.Api;
import jakarta.servlet.http.HttpServletRequest;
import java.time.Instant;
import java.time.temporal.ChronoUnit;
import java.util.HashMap;
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
  private final SamService samService;
  private final PipelinesService pipelinesService;
  private final PipelineRunsService pipelineRunsService;
  private final PipelineInputsOutputsService pipelineInputsOutputsService;
  private final QuotasService quotasService;
  private final DownloadCallCounterService downloadCallCounterService;
  private final IngressConfiguration ingressConfiguration;
  private final PipelineConfigurations pipelinesConfigurations;

  @Autowired
  public PipelineRunsApiController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      JobService jobService,
      SamService samService,
      PipelinesService pipelinesService,
      PipelineRunsService pipelineRunsService,
      PipelineInputsOutputsService pipelineInputsOutputsService,
      QuotasService quotasService,
      DownloadCallCounterService downloadCallCounterService,
      IngressConfiguration ingressConfiguration,
      PipelineConfigurations pipelinesConfigurations) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.jobService = jobService;
    this.samService = samService;
    this.pipelinesService = pipelinesService;
    this.pipelineRunsService = pipelineRunsService;
    this.pipelineInputsOutputsService = pipelineInputsOutputsService;
    this.quotasService = quotasService;
    this.downloadCallCounterService = downloadCallCounterService;
    this.ingressConfiguration = ingressConfiguration;
    this.pipelinesConfigurations = pipelinesConfigurations;
  }

  private static final Logger logger = LoggerFactory.getLogger(PipelineRunsApiController.class);

  private SamUser getAuthenticatedInfo() {
    return samUserFactory.from(request, samConfiguration.baseUri());
  }

  private static final String PIPELINE_RUN_NOT_FOUND_MESSAGE = "Pipeline run %s not found";

  // PipelineRuns
  /**
   * Prepares a pipeline run by validating inputs, generating signed URLs for local file inputs, and
   * storing job metadata in the database. V2 adds support for cloud-based file inputs.
   *
   * @param body the API request body containing inputs for the pipeline
   * @return the prepared pipeline run response, which includes the job ID and signed URLs for
   *     uploading file inputs
   */
  @Override
  public ResponseEntity<ApiPreparePipelineRunResponseV2> preparePipelineRunV2(
      @RequestBody ApiPreparePipelineRunRequestBody body) {
    final SamUser authedUser = getAuthenticatedInfo();
    boolean showHiddenPipelines = samService.isAdmin(authedUser);
    String userId = authedUser.getSubjectId();
    UUID jobId = body.getJobId();
    String pipelineName = body.getPipelineName();
    String description = body.getDescription();
    Boolean useResumableUploads = body.isUseResumableUploads();

    Integer pipelineVersion = body.getPipelineVersion();
    Map<String, Object> userProvidedInputs = body.getPipelineInputs();

    // validate the pipeline name and user-provided inputs
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);
    Pipeline pipeline =
        pipelinesService.getPipeline(validatedPipelineName, pipelineVersion, showHiddenPipelines);

    pipelineInputsOutputsService.validateUserProvidedInputsWithCloud(
        pipeline.getPipelineInputDefinitions(), userProvidedInputs);

    // validate that user has enough quota to run the pipeline
    quotasService.validateUserHasEnoughQuota(userId, validatedPipelineName);

    logger.info(
        "Preparing {} pipeline (version {}) job (id {}) for user {} with validated inputs {}",
        pipelineName,
        pipelineVersion,
        jobId,
        userId,
        userProvidedInputs);

    Map<String, Map<String, String>> fileInputUploadUrls =
        pipelineRunsService.preparePipelineRunV2(
            pipeline, jobId, userId, userProvidedInputs, description, useResumableUploads);

    ApiPreparePipelineRunResponseV2 prepareResponse =
        new ApiPreparePipelineRunResponseV2().jobId(jobId).fileInputUploadUrls(fileInputUploadUrls);

    return new ResponseEntity<>(prepareResponse, HttpStatus.OK);
  }

  /**
   * Prepares a pipeline run by validating inputs, generating signed URLs for local file inputs, and
   * storing job metadata in the database.
   *
   * @param body the API request body containing inputs for the pipeline
   * @return the prepared pipeline run response, which includes the job ID and signed URLs for
   *     uploading file inputs
   * @deprecated use preparePipelineRunV2
   */
  @Deprecated(since = "2.2.0")
  @Override
  public ResponseEntity<ApiPreparePipelineRunResponse> preparePipelineRun(
      @RequestBody ApiPreparePipelineRunRequestBody body) {
    final SamUser authedUser = getAuthenticatedInfo();
    boolean showHiddenPipelines = samService.isAdmin(authedUser);
    String userId = authedUser.getSubjectId();
    UUID jobId = body.getJobId();
    String pipelineName = body.getPipelineName();
    String description = body.getDescription();
    Boolean useResumableUploads = body.isUseResumableUploads();

    Integer pipelineVersion = body.getPipelineVersion();
    Map<String, Object> userProvidedInputs = body.getPipelineInputs();

    // validate the pipeline name and user-provided inputs
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);
    Pipeline pipeline =
        pipelinesService.getPipeline(validatedPipelineName, pipelineVersion, showHiddenPipelines);

    pipelineInputsOutputsService.validateUserProvidedInputs(
        pipeline.getPipelineInputDefinitions(), userProvidedInputs);

    // validate that user has enough quota to run the pipeline
    quotasService.validateUserHasEnoughQuota(userId, validatedPipelineName);

    logger.info(
        "Preparing {} pipeline (version {}) job (id {}) for user {} with validated inputs {}",
        pipelineName,
        pipelineVersion,
        jobId,
        userId,
        userProvidedInputs);

    Map<String, Map<String, String>> fileInputUploadUrls =
        pipelineRunsService.preparePipelineRun(
            pipeline, jobId, userId, userProvidedInputs, description, useResumableUploads);

    ApiPreparePipelineRunResponse prepareResponse =
        new ApiPreparePipelineRunResponse().jobId(jobId).fileInputUploadUrls(fileInputUploadUrls);

    return new ResponseEntity<>(prepareResponse, HttpStatus.OK);
  }

  /**
   * Kicks off the asynchronous process (managed by Stairway) running the specified pipeline job and
   * delivering the outputs to the user
   *
   * @param body the API request body containing the job ID to start
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
        .plus(pipelinesConfigurations.getCommon().getUserDataTtlDays(), ChronoUnit.DAYS)
        .isBefore(Instant.now())) {
      throw new BadRequestException(
          "Pipeline run was prepared more than %s days ago; it cannot be started"
              .formatted(pipelinesConfigurations.getCommon().getUserDataTtlDays()));
    }

    PipelineRun pipelineRunAfterStart =
        pipelineRunsService.startPipelineRun(pipeline, jobId, userId);

    ApiAsyncPipelineRunResponse createdRunResponse =
        pipelineRunToApi(pipelineRunAfterStart, pipeline);

    return new ResponseEntity<>(
        createdRunResponse, getAsyncResponseCode(createdRunResponse.getJobReport()));
  }

  /**
   * Retrieves the result of a pipeline run, including job report, error report (if any), and
   * pipeline run report, including signed urls for outputs if available.
   *
   * @param jobId the ID of the job to retrieve
   * @return the pipeline run result
   * @deprecated use getPipelineRunResultV2
   */
  @Deprecated(since = "2.1.0")
  @Override
  public ResponseEntity<ApiAsyncPipelineRunResponse> getPipelineRunResult(
      @PathVariable("jobId") UUID jobId) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();

    PipelineRun pipelineRun = pipelineRunsService.getPipelineRun(jobId, userId);
    if (pipelineRun == null) {
      throw new NotFoundException(PIPELINE_RUN_NOT_FOUND_MESSAGE.formatted(jobId));
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
   * Retrieves the result of a pipeline run, including job report, error report (if any), and
   * pipeline run report.
   *
   * @param jobId the ID of the job to retrieve
   * @return the pipeline run result
   */
  @Override
  public ResponseEntity<ApiAsyncPipelineRunResponseV2> getPipelineRunResultV2(
      @PathVariable("jobId") UUID jobId) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();

    PipelineRun pipelineRun = pipelineRunsService.getPipelineRun(jobId, userId);
    if (pipelineRun == null) {
      throw new NotFoundException(PIPELINE_RUN_NOT_FOUND_MESSAGE.formatted(jobId));
    }

    if (pipelineRun.getStatus().equals(CommonPipelineRunStatusEnum.PREPARING)) {
      throw new BadRequestException(
          "Pipeline run %s is still preparing; it has to be started before you can query the result"
              .formatted(jobId));
    }

    Pipeline pipeline = pipelinesService.getPipelineById(pipelineRun.getPipelineId());

    ApiAsyncPipelineRunResponseV2 runResponse = pipelineRunToApiV2(pipelineRun, pipeline);

    return new ResponseEntity<>(runResponse, getAsyncResponseCode(runResponse.getJobReport()));
  }

  @Override
  public ResponseEntity<ApiPipelineRunOutputSignedUrlsResponse> getPipelineRunOutputSignedUrls(
      @PathVariable("jobId") UUID jobId) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();

    PipelineRun pipelineRun = pipelineRunsService.getPipelineRun(jobId, userId);
    if (pipelineRun == null) {
      throw new NotFoundException(PIPELINE_RUN_NOT_FOUND_MESSAGE.formatted(jobId));
    }

    if (!pipelineRun.getStatus().equals(CommonPipelineRunStatusEnum.SUCCEEDED)) {
      throw new BadRequestException(
          "Pipeline run %s has state %s; output signed URLs can only be retrieved for complete and successful runs"
              .formatted(jobId, pipelineRun.getStatus()));
    }

    Instant outputExpirationDate = calculateOutputExpirationDate(pipelineRun);
    if (outputExpirationDate.isBefore(Instant.now())) {
      throw new BadRequestException(
          "Outputs for pipeline run %s have expired and are no longer available for download"
              .formatted(jobId));
    }

    ApiPipelineRunOutputSignedUrlsResponse response =
        new ApiPipelineRunOutputSignedUrlsResponse()
            .jobId(jobId)
            .outputSignedUrls(
                pipelineInputsOutputsService.generatePipelineRunOutputSignedUrls(pipelineRun))
            .outputExpirationDate(calculateOutputExpirationDate(pipelineRun).toString());

    downloadCallCounterService.incrementDownloadCallCount(jobId);

    return new ResponseEntity<>(response, HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiGetPipelineRunsResponseV2> getAllPipelineRunsV2(
      Integer pageNumber,
      Integer pageSize,
      String sortProperty,
      String sortDirection,
      String status,
      UUID jobId,
      String pipelineName,
      String description) {
    final SamUser authedUser = getAuthenticatedInfo();
    String userId = authedUser.getSubjectId();
    int maxPageSize = Math.min(pageSize, 100);

    // build filter map from individual parameters
    Map<String, String> filterOptions = new HashMap<>();
    if (status != null) filterOptions.put(PipelineRunFilterSpecification.FILTER_STATUS, status);
    if (jobId != null)
      filterOptions.put(PipelineRunFilterSpecification.FILTER_JOB_ID, jobId.toString());
    if (pipelineName != null)
      filterOptions.put(PipelineRunFilterSpecification.FILTER_PIPELINE_NAME, pipelineName);
    if (description != null)
      filterOptions.put(PipelineRunFilterSpecification.FILTER_DESCRIPTION, description);

    // grab results from current page based on user provided inputs
    // PageRequest is zero-indexed, but for the API we want to be one-indexed for user-friendliness
    Page<PipelineRun> pageResults =
        pipelineRunsService.findPipelineRunsPaginated(
            pageNumber - 1, maxPageSize, sortProperty, sortDirection, userId, filterOptions);

    // convert list of pipelines to map of id to pipeline for all pipelines
    Map<Long, Pipeline> pipelineIdToPipeline =
        pipelinesService.getPipelines(true).stream()
            .collect(Collectors.toMap(Pipeline::getId, p -> p));

    int totalResults = Math.toIntExact(pipelineRunsService.getPipelineRunCount(userId));

    int totalFilteredResults =
        Math.toIntExact(pipelineRunsService.getFilteredPipelineRunCount(userId, filterOptions));

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
                                : null)
                        .outputExpirationDate(
                            pipelineRun.getStatus().isSuccess()
                                ? calculateOutputExpirationDate(pipelineRun).toString()
                                : null))
            .toList();

    ApiGetPipelineRunsResponseV2 apiGetPipelineRunsResponse =
        new ApiGetPipelineRunsResponseV2()
            .results(apiPipelineRuns)
            .totalResults(totalResults)
            .totalFilteredResults(totalFilteredResults);
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
   * @param pipeline the Pipeline associated with the PipelineRun
   * @return ApiAsyncPipelineRunResponse
   */
  private ApiAsyncPipelineRunResponse pipelineRunToApi(PipelineRun pipelineRun, Pipeline pipeline) {
    ApiAsyncPipelineRunResponse response = new ApiAsyncPipelineRunResponse();

    ApiPipelineUserProvidedInputs userProvidedInputs = new ApiPipelineUserProvidedInputs();
    userProvidedInputs.putAll(pipelineInputsOutputsService.retrieveUserProvidedInputs(pipelineRun));

    response.pipelineRunReport(
        new ApiPipelineRunReport()
            .pipelineName(pipeline.getName().getValue())
            .pipelineVersion(pipeline.getVersion())
            .toolVersion(
                pipelineRun
                    .getToolVersion()) // toolVersion comes from pipelineRun, since the pipeline
            // might have been updated since the pipelineRun began
            .userInputs(userProvidedInputs));

    Integer inputSize = pipelineRun.getRawQuotaConsumed();
    if (inputSize != null) {
      String inputSizeUnits = quotasService.getQuotaUnitsForPipeline(pipeline.getName()).getValue();
      response.getPipelineRunReport().inputSize(inputSize).inputSizeUnits(inputSizeUnits);
    }

    // if the pipeline run is successful, return the job report and add outputs to the response
    if (pipelineRun.getStatus().isSuccess()) {
      // calculate the expiration date for the output files
      Instant outputExpirationDate = calculateOutputExpirationDate(pipelineRun);
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
                          ingressConfiguration.getDomainName(), pipelineRun.getJobId(), 1)))
          .pipelineRunReport(
              response
                  .getPipelineRunReport()
                  .outputExpirationDate(outputExpirationDate.toString())
                  .quotaConsumed(pipelineRun.getQuotaConsumed()));

      // return outputs if we have not passed the output expiration date
      if (outputExpirationDate.isAfter(Instant.now())) {
        response
            .getPipelineRunReport()
            .outputs(pipelineInputsOutputsService.generatePipelineRunOutputSignedUrls(pipelineRun));
      }
      return response;

    } else {
      JobApiUtils.AsyncJobResult<String> jobResult =
          jobService.retrieveAsyncJobResult(
              pipelineRun.getJobId(), pipelineRun.getUserId(), String.class, null);
      return response
          .jobReport(jobResult.getJobReport())
          .errorReport(jobResult.getApiErrorReport())
          .pipelineRunReport(
              response.getPipelineRunReport().quotaConsumed(pipelineRun.getQuotaConsumed()));
    }
  }

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
  private ApiAsyncPipelineRunResponseV2 pipelineRunToApiV2(
      PipelineRun pipelineRun, Pipeline pipeline) {
    ApiAsyncPipelineRunResponseV2 response = new ApiAsyncPipelineRunResponseV2();

    ApiPipelineUserProvidedInputs userProvidedInputs = new ApiPipelineUserProvidedInputs();
    userProvidedInputs.putAll(pipelineInputsOutputsService.retrieveUserProvidedInputs(pipelineRun));

    response.pipelineRunReport(
        new ApiPipelineRunReportV2()
            .pipelineName(pipeline.getName().getValue())
            .pipelineVersion(pipeline.getVersion())
            .toolVersion(
                pipelineRun
                    .getToolVersion()) // toolVersion comes from pipelineRun, since the pipeline
            // might have been updated since the pipelineRun began
            .userInputs(userProvidedInputs));

    Integer inputSize = pipelineRun.getRawQuotaConsumed();
    if (inputSize != null) {
      String inputSizeUnits = quotasService.getQuotaUnitsForPipeline(pipeline.getName()).getValue();
      response.getPipelineRunReport().inputSize(inputSize).inputSizeUnits(inputSizeUnits);
    }

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
                          ingressConfiguration.getDomainName(), pipelineRun.getJobId(), 2)))
          .pipelineRunReport(
              response
                  .getPipelineRunReport()
                  .outputs(pipelineInputsOutputsService.getPipelineRunOutputs(pipelineRun))
                  .outputExpirationDate(calculateOutputExpirationDate(pipelineRun).toString())
                  .quotaConsumed(pipelineRun.getQuotaConsumed()));
    } else {
      JobApiUtils.AsyncJobResult<String> jobResult =
          jobService.retrieveAsyncJobResult(
              pipelineRun.getJobId(), pipelineRun.getUserId(), String.class, null);
      return response
          .jobReport(jobResult.getJobReport())
          .errorReport(jobResult.getApiErrorReport())
          .pipelineRunReport(
              response.getPipelineRunReport().quotaConsumed(pipelineRun.getQuotaConsumed()));
    }
  }

  /**
   * Calculate the output expiration date for a pipeline run based on its updated timestamp and the
   * user data TTL configuration.
   *
   * @param pipelineRun the pipeline run
   * @return the output expiration date as an Instant
   */
  private Instant calculateOutputExpirationDate(PipelineRun pipelineRun) {
    return pipelineRun
        .getUpdated()
        .plus(pipelinesConfigurations.getCommon().getUserDataTtlDays(), ChronoUnit.DAYS);
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
