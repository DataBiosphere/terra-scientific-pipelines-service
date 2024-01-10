package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.ApiException;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.common.MetricsUtils;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.exception.InvalidPipelineException;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.generated.api.PipelinesApi;
import bio.terra.pipelines.generated.model.*;
import bio.terra.pipelines.service.ImputationService;
import bio.terra.pipelines.service.PipelinesService;
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
  private final StairwayJobService stairwayJobService;
  private final PipelinesService pipelinesService;
  private final ImputationService imputationService;

  @Autowired
  public PipelinesApiController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      StairwayJobService stairwayJobService,
      PipelinesService pipelinesService,
      ImputationService imputationService) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.pipelinesService = pipelinesService;
    this.stairwayJobService = stairwayJobService;
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
    Pipeline pipelineInfo = pipelinesService.getPipeline(pipelineId);
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

  @Override
  public ResponseEntity<ApiCreateJobResult> createJob(
      @PathVariable("pipelineId") String pipelineId, @RequestBody ApiCreateJobRequestBody body) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    String pipelineVersion = body.getPipelineVersion();
    Object pipelineInputs = body.getPipelineInputs();

    PipelinesEnum validatedPipelineId = validatePipelineId(pipelineId);
    if (validatedPipelineId == null) {
      throw new InvalidPipelineException(String.format("%s is not a valid pipelineId", pipelineId));
    }

    logger.info(
        "Creating {} pipeline job (version {}) for {} user {} with inputs {}",
        pipelineId,
        pipelineVersion,
        userRequest.getEmail(),
        userId,
        pipelineInputs);

    // TODO: make ticket to have the uuid be provided by the caller
    UUID createdJobUuid;
    switch (validatedPipelineId) {
      case IMPUTATION:
        // eventually we'll expand this out to kick off the imputation pipeline flight but for
        // now this is good enough.
        imputationService.queryForWorkspaceApps();

        createdJobUuid =
            imputationService.createImputationJob(userId, pipelineVersion, pipelineInputs);
        break;
      default:
        // this really should never happen since we validate the pipelineId above
        logger.error("Unknown pipeline id {}", pipelineId);
        throw new ApiException("An internal error occurred.");
    }

    if (createdJobUuid == null) {
      logger.error("New {} pipeline job creation failed.", pipelineId);
      throw new ApiException("An internal error occurred.");
    }

    logger.info("Created {} job {}", pipelineId, createdJobUuid);

    ApiJobControl createdJobControl = new ApiJobControl().id(createdJobUuid.toString());
    ApiCreateJobResult createdJobResult = new ApiCreateJobResult().jobControl(createdJobControl);

    MetricsUtils.incrementPipelineRun(pipelineId);

    return new ResponseEntity<>(createdJobResult, HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiGetJobsResponse> getPipelineJobs(
      @PathVariable("pipelineId") String pipelineId, Integer limit, String pageToken) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    validatePipelineId(pipelineId);
    EnumeratedJobs enumeratedJobs =
        stairwayJobService.enumerateJobs(userId, limit, pageToken, pipelineId);

    ApiGetJobsResponse result = JobApiUtils.mapEnumeratedJobsToApi(enumeratedJobs);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  private PipelinesEnum validatePipelineId(String pipelineId) {
    try {
      return PipelinesEnum.valueOf(pipelineId.toUpperCase());
    } catch (IllegalArgumentException e) {
      logger.error("Unknown pipeline id {}", pipelineId);
      return null;
    }
  }
}
