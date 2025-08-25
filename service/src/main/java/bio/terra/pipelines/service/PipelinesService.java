package bio.terra.pipelines.service;

import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.rawls.model.WorkspaceDetails;
import jakarta.validation.constraints.NotNull;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;
import org.springframework.validation.annotation.Validated;

/** The Pipelines Service manages information about the service's available Scientific Pipelines. */
@Component
@Validated
public class PipelinesService {
  private static final Logger logger = LoggerFactory.getLogger(PipelinesService.class);

  private final PipelinesRepository pipelinesRepository;
  private final RawlsService rawlsService;
  private final SamService samService;

  private static final String SEM_VER_REGEX_STRING =
      "^(0|[1-9]\\d*|[a-zA-Z_]*v\\d+)\\.(0|[1-9]\\d*)\\.(0|[1-9]\\d*)$";

  @Autowired
  public PipelinesService(
      PipelinesRepository pipelinesRepository, RawlsService rawlsService, SamService samService) {
    this.pipelinesRepository = pipelinesRepository;
    this.rawlsService = rawlsService;
    this.samService = samService;
  }

  public List<Pipeline> getPipelines() {
    logger.info("Get all Pipelines");
    return pipelinesRepository.findAll();
  }

  public Pipeline getPipeline(PipelinesEnum pipelineName, Integer pipelineVersion) {
    logger.info(
        "Get a specific pipeline for pipelineName {} and version {}",
        pipelineName,
        pipelineVersion);
    if (pipelineVersion == null) {
      return getLatestPipeline(pipelineName);
    }
    Pipeline dbResult = pipelinesRepository.findByNameAndVersion(pipelineName, pipelineVersion);
    if (dbResult == null) {
      throw new NotFoundException(
          "Pipeline not found for pipelineName %s and version %s"
              .formatted(pipelineName, pipelineVersion));
    }
    return dbResult;
  }

  public Pipeline getLatestPipeline(PipelinesEnum pipelineName) {
    logger.info("Get the latest pipeline for pipelineName {}", pipelineName);
    Pipeline dbResult = pipelinesRepository.findFirstByNameOrderByVersionDesc(pipelineName);
    if (dbResult == null) {
      throw new NotFoundException("Pipeline not found for pipelineName %s".formatted(pipelineName));
    }
    return dbResult;
  }

  public Pipeline getPipelineById(Long pipelineId) {
    logger.info("Get a specific pipeline for pipelineId {}", pipelineId);
    Pipeline dbResult = pipelinesRepository.findById(pipelineId).orElse(null);
    if (dbResult == null) {
      throw new NotFoundException("Pipeline not found for pipelineId %s".formatted(pipelineId));
    }
    return dbResult;
  }

  /**
   * This method is meant to only be called by an admin endpoint to update pipeline parameters such
   * as control workspace information and wdl method version.
   *
   * <p>Calls Rawls to fetch control workspace metadata based on the workspaceBillingProject and
   * workspaceName.
   *
   * @param pipelineName - name of pipeline to update
   * @param workspaceBillingProject - workspace billing project to update to
   * @param workspaceName - workspace name to update to
   * @param toolVersion - version of the tool expected to run for corresponding pipeline. must align
   *     with pipeline version
   */
  public Pipeline adminUpdatePipelineWorkspace(
      PipelinesEnum pipelineName,
      Integer pipelineVersion,
      @NotNull String workspaceBillingProject,
      @NotNull String workspaceName,
      @NotNull String toolVersion) {
    WorkspaceDetails workspaceDetails =
        rawlsService.getWorkspaceDetails(
            samService.getTeaspoonsServiceAccountToken(), workspaceBillingProject, workspaceName);
    String workspaceStorageContainerUrl = rawlsService.getWorkspaceBucketName(workspaceDetails);
    String workspaceGoogleProject = rawlsService.getWorkspaceGoogleProject(workspaceDetails);

    Pipeline pipeline = getPipeline(pipelineName, pipelineVersion);
    pipeline.setWorkspaceBillingProject(workspaceBillingProject);
    pipeline.setWorkspaceName(workspaceName);
    pipeline.setWorkspaceStorageContainerName(workspaceStorageContainerUrl);
    pipeline.setWorkspaceGoogleProject(workspaceGoogleProject);

    // ensure toolVersion follows semantic versioning regex (can be preceded by a string ending
    // in v)
    final Pattern pattern = Pattern.compile(SEM_VER_REGEX_STRING);
    final Matcher matcher = pattern.matcher(toolVersion);
    if (!matcher.matches()) {
      throw new ValidationException(
          String.format(
              "toolVersion %s does not follow semantic versioning regex %s",
              toolVersion, SEM_VER_REGEX_STRING));
    }
    pipeline.setToolVersion(toolVersion);

    pipelinesRepository.save(pipeline);
    return pipeline;
  }
}
