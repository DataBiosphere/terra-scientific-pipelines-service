package bio.terra.pipelines.app.controller;

import static bio.terra.pipelines.common.utils.PipelinesEnum.IMPUTATION;

import bio.terra.common.exception.ApiException;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.exception.InvalidPipelineException;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.generated.api.PipelinesApi;
import bio.terra.pipelines.generated.model.*;
import bio.terra.pipelines.service.ImputationService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.stairway.FlightState;
import io.swagger.annotations.Api;
import java.util.List;
import java.util.UUID;
import javax.servlet.http.HttpServletRequest;
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

  @Autowired
  public PipelinesApiController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      JobService jobService,
      PipelinesService pipelinesService,
      ImputationService imputationService) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.pipelinesService = pipelinesService;
    this.jobService = jobService;
    this.imputationService = imputationService;
  }

  private static final Logger logger = LoggerFactory.getLogger(PipelinesApiController.class);

  private SamUser getAuthenticatedInfo() {
    return samUserFactory.from(request, samConfiguration.baseUri());
  }

  // -- Pipelines --

  @Override
  public ResponseEntity<ApiGetPipelinesResult> getPipelines() {
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    ApiGetPipelinesResult result = pipelinesToApi(pipelineList);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiPipeline> getPipeline(@PathVariable("pipelineId") String pipelineId) {
    PipelinesEnum validatedPipelineId = validatePipelineId(pipelineId);
    Pipeline pipelineInfo = pipelinesService.getPipeline(validatedPipelineId);
    ApiPipeline result = pipelineToApi(pipelineInfo);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  static ApiGetPipelinesResult pipelinesToApi(List<Pipeline> pipelineList) {
    ApiGetPipelinesResult apiResult = new ApiGetPipelinesResult();

    for (Pipeline pipeline : pipelineList) {
      apiResult.add(pipelineToApi(pipeline));
    }

    return apiResult;
  }

  static ApiPipeline pipelineToApi(Pipeline pipelineInfo) {
    return new ApiPipeline()
        .pipelineId(pipelineInfo.getPipelineId())
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
   * @param pipelineId the pipeline to run
   * @param body the inputs for the pipeline
   * @return the created job response, which includes a job report containing the job ID,
   *     description, status, status code, submitted timestamp, completed timestamp (if completed),
   *     and result URL. The response also includes an error report if the job failed.
   */
  @Override
  public ResponseEntity<ApiCreateJobResponse> createJob(
      @PathVariable("pipelineId") String pipelineId, @RequestBody ApiCreateJobRequestBody body) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    String description = body.getDescription();
    String pipelineVersion = body.getPipelineVersion();
    Object pipelineInputs = body.getPipelineInputs();

    PipelinesEnum validatedPipelineId = validatePipelineId(pipelineId);

    logger.info(
        "Creating {} pipeline job (version {}) for {} user {} with inputs {}",
        pipelineId,
        pipelineVersion,
        userRequest.getEmail(),
        userId,
        pipelineInputs);

    // TSPS-136 will require that the user provide the job UUID
    UUID createdJobUuid;
    if (validatedPipelineId == IMPUTATION) {
      // eventually we'll expand this out to kick off the imputation pipeline flight but for
      // now this is good enough.
      imputationService.queryForWorkspaceApps();

      createdJobUuid =
          imputationService.createImputationJob(
              userId, description, pipelineVersion, pipelineInputs);
    } else {
      logger.error("Unknown validatedPipelineId {}", validatedPipelineId);
      throw new ApiException("An internal error occurred.");
    }

    logger.info("Created {} job {}", validatedPipelineId.getValue(), createdJobUuid);

    MetricsUtils.incrementPipelineRun(validatedPipelineId);

    FlightState flightState = jobService.retrieveJob(createdJobUuid, userId);
    ApiJobReport jobReport = JobApiUtils.mapFlightStateToApiJobReport(flightState);
    ApiCreateJobResponse createdJobResponse = new ApiCreateJobResponse().jobReport(jobReport);

    return new ResponseEntity<>(createdJobResponse, HttpStatus.valueOf(jobReport.getStatusCode()));
  }

  /** Retrieves job reports for all jobs of the specified pipeline that the user has access to. */
  @Override
  public ResponseEntity<ApiGetJobsResponse> getPipelineJobs(
      @PathVariable("pipelineId") String pipelineId, Integer limit, String pageToken) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    PipelinesEnum validatedPipelineId = validatePipelineId(pipelineId);
    EnumeratedJobs enumeratedJobs =
        jobService.enumerateJobs(userId, limit, pageToken, validatedPipelineId);

    ApiGetJobsResponse result = JobApiUtils.mapEnumeratedJobsToApi(enumeratedJobs);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  /**
   * Validates that the pipelineId is a valid pipelineId and returns the Enum value for the
   * pipelineId
   *
   * <p>Note that in PipelinesServiceTest, we check that all the pipelines in the enum exist in the
   * pipelines table
   *
   * @param pipelineId the pipelineId to validate
   * @return the Enum value for the pipelineId
   * @throws InvalidPipelineException if the pipelineId is not valid
   */
  public PipelinesEnum validatePipelineId(String pipelineId) {
    try {
      return PipelinesEnum.valueOf(pipelineId.toUpperCase());
    } catch (IllegalArgumentException e) {
      logger.error("Unknown pipeline id {}", pipelineId);
      throw new InvalidPipelineException(String.format("%s is not a valid pipelineId", pipelineId));
    }
  }
}
