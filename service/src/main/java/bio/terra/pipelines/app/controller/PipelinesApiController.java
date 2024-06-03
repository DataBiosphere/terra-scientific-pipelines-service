package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.NotFoundException;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.generated.api.PipelinesApi;
import bio.terra.pipelines.generated.model.ApiCreateJobRequestBody;
import bio.terra.pipelines.generated.model.ApiCreatePipelineRunResponse;
import bio.terra.pipelines.generated.model.ApiGetPipelinesResult;
import bio.terra.pipelines.generated.model.ApiJobReport;
import bio.terra.pipelines.generated.model.ApiPipeline;
import bio.terra.pipelines.generated.model.ApiPipelineUserProvidedInputDefinition;
import bio.terra.pipelines.generated.model.ApiPipelineUserProvidedInputDefinitions;
import bio.terra.pipelines.generated.model.ApiPipelineWithDetails;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.PipelinesService;
import io.swagger.annotations.Api;
import jakarta.servlet.http.HttpServletRequest;
import java.nio.file.Path;
import java.util.HashMap;
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

/** Pipelines controller */
@Controller
@Api(tags = {"pipelines"})
public class PipelinesApiController implements PipelinesApi {
  private final SamConfiguration samConfiguration;
  private final SamUserFactory samUserFactory;
  private final HttpServletRequest request;
  private final JobService jobService;
  private final PipelinesService pipelinesService;
  private final PipelineRunsService pipelineRunsService;
  private final IngressConfiguration ingressConfiguration;

  @Autowired
  public PipelinesApiController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      JobService jobService,
      PipelinesService pipelinesService,
      PipelineRunsService pipelineRunsService,
      IngressConfiguration ingressConfiguration) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.pipelinesService = pipelinesService;
    this.pipelineRunsService = pipelineRunsService;
    this.jobService = jobService;
    this.ingressConfiguration = ingressConfiguration;
  }

  private static final Logger logger = LoggerFactory.getLogger(PipelinesApiController.class);

  private SamUser getAuthenticatedInfo() {
    return samUserFactory.from(request, samConfiguration.baseUri());
  }

  // -- Pipelines --

  @Override
  public ResponseEntity<ApiGetPipelinesResult> getPipelines() {
    getAuthenticatedInfo();
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    ApiGetPipelinesResult result = pipelinesToApi(pipelineList);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiPipelineWithDetails> getPipelineDetails(
      @PathVariable("pipelineName") String pipelineName) {
    getAuthenticatedInfo();
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);
    Pipeline pipelineInfo = pipelinesService.getPipeline(validatedPipelineName);
    ApiPipelineWithDetails result = pipelineWithDetailsToApi(pipelineInfo);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  static ApiGetPipelinesResult pipelinesToApi(List<Pipeline> pipelineList) {
    ApiGetPipelinesResult apiResult = new ApiGetPipelinesResult();

    for (Pipeline pipeline : pipelineList) {
      apiResult.add(pipelineToApi(pipeline));
    }

    return apiResult;
  }

  static ApiPipelineWithDetails pipelineWithDetailsToApi(Pipeline pipelineInfo) {
    ApiPipelineUserProvidedInputDefinitions inputs = new ApiPipelineUserProvidedInputDefinitions();
    inputs.addAll(
        pipelineInfo.getPipelineInputDefinitions().stream()
            .filter(PipelineInputDefinition::getUserProvided)
            .map(
                input ->
                    new ApiPipelineUserProvidedInputDefinition()
                        .name(input.getName())
                        .type(input.getType().toString())
                        .isRequired(input.getIsRequired()))
            .toList());
    return new ApiPipelineWithDetails()
        .pipelineName(pipelineInfo.getName())
        .displayName(pipelineInfo.getDisplayName())
        .description(pipelineInfo.getDescription())
        .type(pipelineInfo.getPipelineType())
        .inputs(inputs);
  }

  static ApiPipeline pipelineToApi(Pipeline pipelineInfo) {
    return new ApiPipeline()
        .pipelineName(pipelineInfo.getName())
        .displayName(pipelineInfo.getDisplayName())
        .description(pipelineInfo.getDescription());
  }

  // Pipelines runs

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
  public ResponseEntity<ApiCreatePipelineRunResponse> createPipelineRun(
      @PathVariable("pipelineName") String pipelineName,
      @RequestBody ApiCreateJobRequestBody body) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    UUID jobId = body.getJobControl().getId();

    String description = body.getDescription();
    String pipelineVersion = body.getPipelineVersion();
    Map<String, Object> userProvidedInputs = new HashMap<>(body.getPipelineInputs());

    // validate the pipeline name and user-provided inputs
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);
    Pipeline pipeline = pipelinesService.getPipeline(validatedPipelineName);

    pipelinesService.validateUserProvidedInputs(
        pipeline.getPipelineInputDefinitions(), userProvidedInputs);

    logger.info(
        "Creating {} pipeline (version {}) job (id {}) for user {} with validated inputs {}",
        pipelineName,
        pipelineVersion,
        jobId,
        userId,
        userProvidedInputs);

    String resultPath = getAsyncResultEndpoint(ingressConfiguration, request, jobId);

    PipelineRun pipelineRun =
        pipelineRunsService.createPipelineRun(
            pipeline, jobId, userId, description, userProvidedInputs, resultPath);

    ApiCreatePipelineRunResponse createdRunResponse = pipelineRunToApi(pipelineRun);

    return new ResponseEntity<>(
        createdRunResponse, getAsyncResponseCode(createdRunResponse.getJobReport()));
  }

  @Override
  public ResponseEntity<ApiCreatePipelineRunResponse> getPipelineRunResult(
      @PathVariable("pipelineName") String pipelineName, @PathVariable("jobId") UUID jobId) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();

    PipelineApiUtils.validatePipelineName(pipelineName, logger);

    PipelineRun pipelineRun = pipelineRunsService.getPipelineRun(jobId, userId);
    if (pipelineRun == null) {
      throw new NotFoundException("Pipeline run %s not found".formatted(jobId));
    }

    ApiCreatePipelineRunResponse runResponse = pipelineRunToApi(pipelineRun);

    return new ResponseEntity<>(runResponse, getAsyncResponseCode(runResponse.getJobReport()));
  }

  /**
   * Converts a PipelineRun to an ApiCreatePipelineRunResponse. If the PipelineRun has completed
   * successfully (isSuccess is true), we know there are no errors to retrieve and all the
   * information needed to return to the user is available in the pipeline_runs table.
   *
   * <p>If the PipelineRun is not marked as a success, we retrieve the running and/or error
   * information from Stairway.
   *
   * @param pipelineRun the PipelineRun to convert
   * @return ApiCreatePipelineRunResponse
   */
  private ApiCreatePipelineRunResponse pipelineRunToApi(PipelineRun pipelineRun) {
    if (Boolean.TRUE.equals(pipelineRun.getIsSuccess())) {
      return new ApiCreatePipelineRunResponse()
          .jobReport(
              new ApiJobReport()
                  .id(pipelineRun.getJobId().toString())
                  .description(pipelineRun.getDescription())
                  .status(ApiJobReport.StatusEnum.SUCCEEDED)
                  .statusCode(HttpStatus.OK.value())
                  .submitted(pipelineRun.getCreated().toString())
                  .completed(pipelineRun.getUpdated().toString())
                  .resultURL(pipelineRun.getResultUrl()))
          .pipelineOutput(pipelineRun.getOutput());
    } else {
      JobApiUtils.AsyncJobResult<String> jobResult =
          jobService.retrieveAsyncJobResult(
              pipelineRun.getJobId(), pipelineRun.getUserId(), String.class, null);

      return new ApiCreatePipelineRunResponse()
          // TODO - the following will populate the submitted and completed fields with timestamps
          // from stairway, not from the pipeline_runs table. is it worth overwriting those fields?
          .jobReport(jobResult.getJobReport())
          .errorReport(jobResult.getApiErrorReport());
    }
  }

  /**
   * Returns the result endpoint corresponding to an async request. The endpoint is used to build an
   * ApiJobReport. This method retrieves the protocol and domain name from the request and generates
   * a result endpoint with the form: {protocol}{domainName}/{servletpath}/result/{jobId} relative
   * to the async endpoint.
   *
   * @param ingressConfiguration configuration specifying the ingress domain name
   * @param jobId identifier for the job
   * @return a string with the result endpoint URL
   */
  public static String getAsyncResultEndpoint(
      IngressConfiguration ingressConfiguration, HttpServletRequest request, UUID jobId) {
    String endpointPath = "%s/result/%s".formatted(request.getServletPath(), jobId);

    // This is a little hacky, but GCP rejects non-https traffic and a local server does not
    // support it.
    String domainName = ingressConfiguration.getDomainName();
    String protocol = domainName.startsWith("localhost") ? "http://" : "https://";

    return protocol + Path.of(domainName, endpointPath);
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
