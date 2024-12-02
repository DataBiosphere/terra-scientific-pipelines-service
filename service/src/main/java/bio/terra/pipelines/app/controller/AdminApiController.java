package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.generated.api.AdminApi;
import bio.terra.pipelines.generated.model.*;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.service.QuotasService;
import io.swagger.annotations.Api;
import jakarta.servlet.http.HttpServletRequest;
import java.util.Optional;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;

/** Admin controller */
@Controller
@Api(tags = {"admin"})
public class AdminApiController implements AdminApi {
  private final SamConfiguration samConfiguration;
  private final SamUserFactory samUserFactory;
  private final HttpServletRequest request;
  private final PipelinesService pipelinesService;
  private final SamService samService;
  private final QuotasService quotasService;

  @Autowired
  public AdminApiController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      PipelinesService pipelinesService,
      SamService samService,
      QuotasService quotasService) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.pipelinesService = pipelinesService;
    this.samService = samService;
    this.quotasService = quotasService;
  }

  private static final Logger logger = LoggerFactory.getLogger(AdminApiController.class);

  private SamUser getAuthenticatedInfo() {
    return samUserFactory.from(request, samConfiguration.baseUri());
  }

  @Override
  public ResponseEntity<ApiAdminPipeline> getPipeline(
      String pipelineName, Integer pipelineVersion) {
    final SamUser authedUser = getAuthenticatedInfo();
    samService.checkAdminAuthz(authedUser);
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);
    Pipeline pipeline = pipelinesService.getPipeline(validatedPipelineName, pipelineVersion);
    return new ResponseEntity<>(pipelineToApiAdminPipeline(pipeline), HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiAdminQuota> getQuotaForPipelineAndUser(
      String pipelineName, String userId) {
    final SamUser authedUser = getAuthenticatedInfo();
    samService.checkAdminAuthz(authedUser);
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);

    // Check if a row exists for this user + pipeline. Throw an error if it doesn't
    // A row should exist if a user has run into a quota issue before.
    Optional<UserQuota> userQuota =
        quotasService.getQuotaForUserAndPipeline(userId, validatedPipelineName);
    if (userQuota.isEmpty()) {
      throw new BadRequestException(
          String.format(
              "User quota not found for user %s and pipeline %s",
              userId, validatedPipelineName.getValue()));
    }

    return new ResponseEntity<>(userQuotaToApiAdminQuota(userQuota.get()), HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiAdminPipeline> updatePipeline(
      String pipelineName, Integer pipelineVersion, ApiUpdatePipelineRequestBody body) {
    final SamUser authedUser = getAuthenticatedInfo();
    samService.checkAdminAuthz(authedUser);
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);
    String workspaceBillingProject = body.getWorkspaceBillingProject();
    String workspaceName = body.getWorkspaceName();
    String wdlMethodVersion = body.getWdlMethodVersion();
    Pipeline updatedPipeline =
        pipelinesService.adminUpdatePipelineWorkspace(
            validatedPipelineName,
            pipelineVersion,
            workspaceBillingProject,
            workspaceName,
            wdlMethodVersion);
    return new ResponseEntity<>(pipelineToApiAdminPipeline(updatedPipeline), HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiAdminQuota> updateQuotaLimitForPipelineAndUser(
      String pipelineName, String userId, ApiUpdateQuotaLimitRequestBody body) {
    final SamUser authedUser = getAuthenticatedInfo();
    samService.checkAdminAuthz(authedUser);
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);

    // Check if a row exists for this user + pipeline. Throw an error if it doesn't
    // A row should exist if a user has run into a quota issue before.
    Optional<UserQuota> userQuota =
        quotasService.getQuotaForUserAndPipeline(userId, validatedPipelineName);
    if (userQuota.isEmpty()) {
      throw new BadRequestException(
          String.format(
              "User quota not found for user %s and pipeline %s",
              userId, validatedPipelineName.getValue()));
    }
    int newQuotaLimit = body.getQuotaLimit();
    UserQuota updatedUserQuota =
        quotasService.adminUpdateQuotaLimit(userQuota.get(), newQuotaLimit);
    return new ResponseEntity<>(userQuotaToApiAdminQuota(updatedUserQuota), HttpStatus.OK);
  }

  public ApiAdminPipeline pipelineToApiAdminPipeline(Pipeline pipeline) {
    return new ApiAdminPipeline()
        .pipelineName(pipeline.getName().getValue())
        .displayName(pipeline.getDisplayName())
        .description(pipeline.getDescription())
        .workspaceBillingProject(pipeline.getWorkspaceBillingProject())
        .workspaceName(pipeline.getWorkspaceName())
        .workspaceStorageContainerName(pipeline.getWorkspaceStorageContainerName())
        .workspaceGoogleProject(pipeline.getWorkspaceGoogleProject())
        .wdlMethodVersion(pipeline.getWdlMethodVersion());
  }

  public ApiAdminQuota userQuotaToApiAdminQuota(UserQuota userQuota) {
    return new ApiAdminQuota()
        .userId(userQuota.getUserId())
        .pipelineName(userQuota.getPipelineName().getValue())
        .quotaLimit(userQuota.getQuota())
        .quotaConsumed(userQuota.getQuotaConsumed());
  }
}
