package bio.terra.pipelines.app.controller;

import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.generated.api.AdminApi;
import bio.terra.pipelines.generated.model.*;
import bio.terra.pipelines.service.PipelinesService;
import io.swagger.annotations.Api;
import jakarta.servlet.http.HttpServletRequest;
import java.util.UUID;
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

  @Autowired
  public AdminApiController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      PipelinesService pipelinesService,
      SamService samService) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.pipelinesService = pipelinesService;
    this.samService = samService;
  }

  private static final Logger logger = LoggerFactory.getLogger(AdminApiController.class);

  private SamUser getAuthenticatedInfo() {
    return samUserFactory.from(request, samConfiguration.baseUri());
  }

  @Override
  public ResponseEntity<ApiAdminPipeline> updatePipelineWorkspaceId(
      String pipelineName, UUID workspaceId) {
    final SamUser authedUser = getAuthenticatedInfo();
    samService.checkAdminAuthz(authedUser);
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);
    Pipeline updatedPipeline =
        pipelinesService.updatePipelineWorkspaceId(validatedPipelineName, workspaceId);
    return new ResponseEntity<>(pipelineToApiAdminPipeline(updatedPipeline), HttpStatus.OK);
  }

  public ApiAdminPipeline pipelineToApiAdminPipeline(Pipeline pipeline) {
    return new ApiAdminPipeline()
        .pipelineName(pipeline.getName())
        .displayName(pipeline.getDisplayName())
        .description(pipeline.getDescription())
        .workspaceId(pipeline.getWorkspaceId());
  }
}
