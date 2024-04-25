package bio.terra.pipelines.app.controller;

import static bio.terra.pipelines.app.controller.JobApiUtils.mapEnumeratedJobsToApi;
import static bio.terra.pipelines.app.controller.JobApiUtils.mapFlightStateToApiJobReport;
import static bio.terra.pipelines.common.utils.PipelinesEnum.IMPUTATION_BEAGLE;

import bio.terra.common.exception.ApiException;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.generated.api.PipelinesApi;
import bio.terra.pipelines.generated.model.*;
import bio.terra.pipelines.service.ImputationService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.stairway.FlightState;
import io.swagger.annotations.Api;
import jakarta.servlet.http.HttpServletRequest;
import java.nio.file.Path;
import java.util.List;
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
  private final ImputationService imputationService;
  private final IngressConfiguration ingressConfiguration;

  @Autowired
  public PipelinesApiController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      JobService jobService,
      PipelinesService pipelinesService,
      ImputationService imputationService,
      IngressConfiguration ingressConfiguration) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.pipelinesService = pipelinesService;
    this.jobService = jobService;
    this.imputationService = imputationService;
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
            .map(
                input ->
                    new ApiPipelineUserProvidedInputDefinition()
                        .name(input.getName())
                        .type(input.getType())
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

  // Pipelines jobs

  /**
   * Kicks off the asynchronous process (managed by Stairway) of gathering user-provided inputs,
   * running the specified pipeline, and delivering the outputs to the user.
   *
   * <p>For now, the job will be created with a random UUID. In the future (TSPS-136), we will
   * require the user to provide a job UUID.
   *
   * @param pipelineName the pipeline to run
   * @param body the inputs for the pipeline
   * @return the created job response, which includes a job report containing the job ID,
   *     description, status, status code, submitted timestamp, completed timestamp (if completed),
   *     and result URL. The response also includes an error report if the job failed.
   */
  @Override
  public ResponseEntity<ApiCreateJobResponse> createJob(
      @PathVariable("pipelineName") String pipelineName,
      @RequestBody ApiCreateJobRequestBody body) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    UUID jobId = body.getJobControl().getId();

    String description = body.getDescription();
    String pipelineVersion = body.getPipelineVersion();
    Object pipelineInputs = body.getPipelineInputs();

    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);

    pipelinesService.validateInputs(validatedPipelineName, pipelineInputs);

    logger.info(
        "Creating {} pipeline (version {}) job (id {}) for user {} with inputs {}",
        pipelineName,
        pipelineVersion,
        jobId,
        userId,
        pipelineInputs);

    String resultPath = getAsyncResultEndpoint(ingressConfiguration, request, jobId);

    if (validatedPipelineName == IMPUTATION_BEAGLE) {
      Pipeline pipeline = pipelinesService.getPipeline(IMPUTATION_BEAGLE);

      imputationService.createImputationJob(
          jobId, userId, description, pipeline, pipelineInputs, resultPath);
    } else {
      logger.error("Unknown validatedPipelineName {}", validatedPipelineName);
      throw new ApiException("An internal error occurred.");
    }

    logger.info("Created {} job {}", validatedPipelineName.getValue(), jobId);

    FlightState flightState = jobService.retrieveJob(jobId, userId, validatedPipelineName);
    ApiJobReport jobReport = mapFlightStateToApiJobReport(flightState);
    ApiCreateJobResponse createdJobResponse = new ApiCreateJobResponse().jobReport(jobReport);

    return new ResponseEntity<>(createdJobResponse, HttpStatus.valueOf(jobReport.getStatusCode()));
  }

  /** Retrieves job reports for all jobs of the specified pipeline that the user has access to. */
  @Override
  public ResponseEntity<ApiGetJobsResponse> getPipelineJobs(
      @PathVariable("pipelineName") String pipelineName, Integer limit, String pageToken) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);
    EnumeratedJobs enumeratedJobs =
        jobService.enumerateJobs(userId, limit, pageToken, validatedPipelineName);

    ApiGetJobsResponse result = mapEnumeratedJobsToApi(enumeratedJobs);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiCreateJobResponse> getPipelineJobResult(
      @PathVariable("pipelineName") String pipelineName, @PathVariable("jobId") UUID jobId) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    PipelinesEnum validatedPipelineId = PipelineApiUtils.validatePipelineName(pipelineName, logger);

    JobApiUtils.AsyncJobResult<String> jobResult =
        jobService.retrieveAsyncJobResult(jobId, userId, validatedPipelineId, String.class, null);

    ApiCreateJobResponse response =
        new ApiCreateJobResponse()
            .jobReport(jobResult.getJobReport())
            .errorReport(jobResult.getApiErrorReport())
            .pipelineOutput(jobResult.getResult()); // this is null unless the job succeeded

    return new ResponseEntity<>(response, getAsyncResponseCode(response.getJobReport()));
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
    String endpointPath = String.format("%s/result/%s", request.getServletPath(), jobId);

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
