package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.params.provider.Arguments.arguments;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.rawls.model.WorkspaceDetails;
import jakarta.validation.ConstraintViolationException;
import java.util.List;
import java.util.stream.Stream;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import org.mockito.InjectMocks;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.bean.override.mockito.MockitoBean;

class PipelinesServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks PipelinesService pipelinesService;
  @Autowired PipelinesRepository pipelinesRepository;
  @MockitoBean SamService samService;
  @MockitoBean RawlsService rawlsService;
  Integer currentPipelineVersion = 1;

  @Test
  void getCorrectNumberOfPipelines() {
    // migrations insert one pipeline (imputation) so make sure we find it
    List<Pipeline> pipelineList = pipelinesService.getPipelines(false);
    assertEquals(1, pipelineList.size());
    String workspaceBillingProject = "testTerraProject";
    String workspaceName = "testTerraWorkspaceName";
    String workspaceStorageContainerName = "testWorkspaceStorageContainerUrl";
    String workspaceGoogleProject = "testWorkspaceGoogleProject";

    // save a new version of the same pipeline
    pipelinesRepository.save(
        new Pipeline(
            PipelinesEnum.ARRAY_IMPUTATION,
            currentPipelineVersion + 1,
            false,
            "pipelineDisplayName",
            "description",
            "pipelineType",
            "wdlUrl",
            "toolName",
            "1.2.1",
            workspaceBillingProject,
            workspaceName,
            workspaceStorageContainerName,
            workspaceGoogleProject,
            null,
            null));

    pipelineList = pipelinesService.getPipelines(false);
    assertEquals(2, pipelineList.size());
    Pipeline savedPipeline = pipelineList.get(1);
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, savedPipeline.getName());
    assertEquals(currentPipelineVersion + 1, savedPipeline.getVersion());
    assertEquals("pipelineDisplayName", savedPipeline.getDisplayName());
    assertEquals("description", savedPipeline.getDescription());
    assertEquals("pipelineType", savedPipeline.getPipelineType());
    assertEquals("wdlUrl", savedPipeline.getWdlUrl());
    assertEquals("toolName", savedPipeline.getToolName());
    assertEquals("1.2.1", savedPipeline.getToolVersion());
    assertEquals(workspaceBillingProject, savedPipeline.getWorkspaceBillingProject());
    assertEquals(workspaceName, savedPipeline.getWorkspaceName());
    assertEquals(workspaceStorageContainerName, savedPipeline.getWorkspaceStorageContainerName());
    assertEquals(workspaceGoogleProject, savedPipeline.getWorkspaceGoogleProject());

    // save a hidden pipeline
    pipelinesRepository.save(
        new Pipeline(
            PipelinesEnum.ARRAY_IMPUTATION,
            currentPipelineVersion + 2,
            true,
            "pipelineDisplayName",
            "description",
            "pipelineType",
            "wdlUrl",
            "toolName",
            "1.2.1",
            workspaceBillingProject,
            workspaceName,
            workspaceStorageContainerName,
            workspaceGoogleProject,
            null,
            null));

    // make sure hidden pipeline is not returned when showHidden is false
    pipelineList = pipelinesService.getPipelines(false);
    assertEquals(2, pipelineList.size());

    pipelineList = pipelinesService.getPipelines(true);
    assertEquals(3, pipelineList.size());
  }

  @Test
  void getPipelineById() {
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;
    Pipeline p = pipelinesService.getPipeline(imputationPipeline, null, false);
    Long pipelineId = p.getId();

    Pipeline pById = pipelinesService.getPipelineById(pipelineId);
    assertEquals(p, pById);

    assertThrows(NotFoundException.class, () -> pipelinesService.getPipelineById(999L));
  }

  @Test
  void getHiddenPipeline() {
    // save a hidden pipeline
    pipelinesRepository.save(
        new Pipeline(
            PipelinesEnum.ARRAY_IMPUTATION,
            currentPipelineVersion + 1,
            true,
            "pipelineDisplayName",
            "description",
            "pipelineType",
            "wdlUrl",
            "toolName",
            "1.2.1",
            "meh",
            "doesnt",
            "matter",
            "probalby",
            null,
            null));

    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;

    // should not be able to grab hidden pipeline when showHidden is false
    assertThrows(
        NotFoundException.class,
        () -> pipelinesService.getPipeline(imputationPipeline, currentPipelineVersion + 1, false));

    Pipeline p = pipelinesService.getPipeline(imputationPipeline, currentPipelineVersion + 1, true);
    assertEquals(currentPipelineVersion + 1, p.getVersion());

    // should grab non-hidden pipeline when pipeline version is not provided even if showHidden is
    // true
    p = pipelinesService.getPipeline(imputationPipeline, null, true);
    assertEquals(currentPipelineVersion, p.getVersion());
  }

  @Test
  void getPipelineByNameAndVersion() {
    // save a new version of the same pipeline that exists in the table
    pipelinesRepository.save(
        new Pipeline(
            PipelinesEnum.ARRAY_IMPUTATION,
            currentPipelineVersion + 1,
            false,
            "pipelineDisplayName",
            "description",
            "pipelineType",
            "wdlUrl",
            "toolName",
            "1.2.1",
            null,
            null,
            null,
            null,
            null,
            null));
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;
    // this should return the highest version of the pipeline
    Pipeline nullVersionPipeline = pipelinesService.getPipeline(imputationPipeline, null, false);
    assertEquals(currentPipelineVersion + 1, nullVersionPipeline.getVersion());

    // this should return the specific version of the pipeline that exists
    Pipeline specificVersionPipeline =
        pipelinesService.getPipeline(imputationPipeline, currentPipelineVersion, false);
    assertEquals(currentPipelineVersion, specificVersionPipeline.getVersion());

    // if asking for unknown version pipeline combo, throw exception
    assertThrows(
        NotFoundException.class,
        () -> pipelinesService.getPipeline(imputationPipeline, 999, false));
  }

  @Test
  void getLatestPipeline() {
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;
    // this should return the highest version of the pipeline
    Pipeline getLatestPipeline = pipelinesService.getLatestPipeline(imputationPipeline);
    assertEquals(currentPipelineVersion, getLatestPipeline.getVersion());

    // save a new version of the same pipeline that exists in the table
    pipelinesRepository.save(
        new Pipeline(
            PipelinesEnum.ARRAY_IMPUTATION,
            100,
            false,
            "pipelineDisplayName",
            "description",
            "pipelineType",
            "wdlUrl",
            "toolName",
            "1.2.1",
            null,
            null,
            null,
            null,
            null,
            null));

    // this should return the new highest version of the pipeline
    getLatestPipeline = pipelinesService.getLatestPipeline(imputationPipeline);
    assertEquals(100, getLatestPipeline.getVersion());
  }

  @Test
  void testPipelineToString() {
    // test .ToString() method on Pipeline Entity
    List<Pipeline> pipelineList = pipelinesService.getPipelines(false);
    assertEquals(1, pipelineList.size());
    for (Pipeline p : pipelineList) {
      assertEquals(
          String.format(
              "Pipeline[pipelineName=%s, version=%s, hidden=%s, displayName=%s, description=%s, pipelineType=%s, wdlUrl=%s, toolName=%s, toolVersion=%s, workspaceBillingProject=%s, workspaceName=%s, workspaceStorageContainerName=%s, workspaceGoogleProject=%s]",
              p.getName(),
              p.getVersion(),
              p.isHidden(),
              p.getDisplayName(),
              p.getDescription(),
              p.getPipelineType(),
              p.getWdlUrl(),
              p.getToolName(),
              p.getToolVersion(),
              p.getWorkspaceBillingProject(),
              p.getWorkspaceName(),
              p.getWorkspaceStorageContainerName(),
              p.getWorkspaceGoogleProject()),
          p.toString());
    }
  }

  @Test
  void testPipelineHashCode() {
    List<Pipeline> pipelineList = pipelinesService.getPipelines(false);
    assertEquals(1, pipelineList.size());
    for (Pipeline p : pipelineList) {
      // 17 and 31 are hardcoded in this hashCode method of this class
      assertEquals(
          new HashCodeBuilder(17, 31)
              .append(p.getId())
              .append(p.getName())
              .append(p.getVersion())
              .append(p.isHidden())
              .append(p.getDisplayName())
              .append(p.getDescription())
              .append(p.getPipelineType())
              .append(p.getWdlUrl())
              .append(p.getToolName())
              .append(p.getToolVersion())
              .append(p.getWorkspaceBillingProject())
              .append(p.getWorkspaceName())
              .append(p.getWorkspaceStorageContainerName())
              .append(p.getWorkspaceGoogleProject())
              .toHashCode(),
          p.hashCode());
    }
  }

  @Test
  void adminUpdatePipelineWorkspace() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.ARRAY_IMPUTATION;
    Pipeline p = pipelinesService.getPipeline(pipelinesEnum, null, true);
    String newWorkspaceBillingProject = "newTestTerraProject";
    String newWorkspaceName = "newTestTerraWorkspaceName";
    String newToolVersion = "0.13.1";
    boolean newHidden = true;

    String newWorkspaceStorageContainerName = "newTestWorkspaceStorageContainerUrl";
    String newWorkspaceGoogleProject = "newTestWorkspaceGoogleProject";
    WorkspaceDetails workspaceDetails =
        new WorkspaceDetails()
            .bucketName(newWorkspaceStorageContainerName)
            .googleProject(newWorkspaceGoogleProject);

    // make sure the current pipeline does not have the workspace info we're trying to update with
    assertNotEquals(newWorkspaceBillingProject, p.getWorkspaceBillingProject());
    assertNotEquals(newWorkspaceName, p.getWorkspaceName());
    assertNotEquals(newWorkspaceStorageContainerName, p.getWorkspaceStorageContainerName());
    assertNotEquals(newWorkspaceGoogleProject, p.getWorkspaceGoogleProject());
    assertNotEquals(newToolVersion, p.getToolVersion());
    assertNotEquals(newHidden, p.isHidden());

    // mocks
    when(samService.getTeaspoonsServiceAccountToken()).thenReturn("fakeToken");
    when(rawlsService.getWorkspaceDetails(
            "fakeToken", newWorkspaceBillingProject, newWorkspaceName))
        .thenReturn(workspaceDetails);
    when(rawlsService.getWorkspaceBucketName(workspaceDetails))
        .thenReturn(newWorkspaceStorageContainerName);
    when(rawlsService.getWorkspaceGoogleProject(workspaceDetails))
        .thenReturn(newWorkspaceGoogleProject);

    // update pipeline values
    pipelinesService.adminUpdatePipelineWorkspace(
        pipelinesEnum,
        p.getVersion(),
        true,
        newWorkspaceBillingProject,
        newWorkspaceName,
        newToolVersion);
    p = pipelinesService.getPipeline(pipelinesEnum, p.getVersion(), true);

    // assert the workspace info has been updated
    assertEquals(newWorkspaceBillingProject, p.getWorkspaceBillingProject());
    assertEquals(newWorkspaceName, p.getWorkspaceName());
    assertEquals(newToolVersion, p.getToolVersion());

    // assert the fetched info from Rawls has been written
    assertEquals(newWorkspaceStorageContainerName, p.getWorkspaceStorageContainerName());
    assertEquals(newWorkspaceGoogleProject, p.getWorkspaceGoogleProject());

    // assert pipeline visibility has been updated
    assertTrue(p.isHidden());
    assertThrows(
        NotFoundException.class, () -> pipelinesService.getPipeline(pipelinesEnum, null, false));

    // update pipeline with null isHidden value - should not change visibility
    pipelinesService.adminUpdatePipelineWorkspace(
        pipelinesEnum,
        p.getVersion(),
        null,
        newWorkspaceBillingProject,
        newWorkspaceName,
        newToolVersion);
    p = pipelinesService.getPipeline(pipelinesEnum, p.getVersion(), true);
    assertTrue(p.isHidden());

    pipelinesService.adminUpdatePipelineWorkspace(
        pipelinesEnum,
        p.getVersion(),
        false,
        newWorkspaceBillingProject,
        newWorkspaceName,
        newToolVersion);
    p = pipelinesService.getPipeline(pipelinesEnum, p.getVersion(), true);
    // assert pipeline visibility has not been updated
    assertFalse(p.isHidden());
  }

  private static Stream<Arguments> badToolVersions() {
    return Stream.of(
        arguments("blah.3.2"),
        arguments("0.2"),
        arguments("0.bhmm.2"),
        arguments("0.4.ok"),
        arguments("0.v3.4"),
        arguments("strv.0.1.2ing"),
        arguments("ImputationBeagle-development_v0.0.0"), // dashes not allowed
        arguments(
            "ImputationBeagle_development_0.0.0"), // missing a "v" before the semantic version
        arguments("hiiv.1.4"));
  }

  @ParameterizedTest
  @MethodSource("badToolVersions")
  void adminUpdatePipelineWorkspaceBadToolVersion(String badToolVersion) {
    PipelinesEnum pipelinesEnum = PipelinesEnum.ARRAY_IMPUTATION;
    Pipeline p = pipelinesService.getPipeline(pipelinesEnum, null, false);
    String newWorkspaceBillingProject = "newTestTerraProject";
    String newWorkspaceName = "newTestTerraWorkspaceName";
    int version = p.getVersion();
    assertThrows(
        ValidationException.class,
        () ->
            pipelinesService.adminUpdatePipelineWorkspace(
                pipelinesEnum,
                version,
                null,
                newWorkspaceBillingProject,
                newWorkspaceName,
                badToolVersion));
  }

  private static Stream<Arguments> goodToolVersions() {
    return Stream.of(
        arguments("0.13.1"),
        arguments("ImputationBeagle_development_v0.0.0"),
        arguments("stringwithvinthemiddlev0.0.0"),
        arguments("v0.1.32"),
        arguments("stringv0.1.32"),
        arguments("1.13.1"));
  }

  @ParameterizedTest
  @MethodSource("goodToolVersions")
  void adminUpdatePipelineWorkspaceGoodToolVersion(String goodToolVersion) {
    PipelinesEnum pipelinesEnum = PipelinesEnum.ARRAY_IMPUTATION;
    Pipeline p = pipelinesService.getPipeline(pipelinesEnum, null, false);
    String newWorkspaceBillingProject = "newTestTerraProject";
    String newWorkspaceName = "newTestTerraWorkspaceName";
    assertDoesNotThrow(
        () ->
            pipelinesService.adminUpdatePipelineWorkspace(
                pipelinesEnum,
                p.getVersion(),
                null,
                newWorkspaceBillingProject,
                newWorkspaceName,
                goodToolVersion));
  }

  @Test
  void adminUpdatePipelineWorkspaceNullValueThrows() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.ARRAY_IMPUTATION;
    Pipeline p = pipelinesService.getPipeline(pipelinesEnum, null, false);

    String newWorkspaceName = "newTestTerraWorkspaceName";
    int version = p.getVersion();

    // make sure the current pipeline does not have the workspace info we're trying to update with
    assertNotEquals(newWorkspaceName, p.getWorkspaceName());

    // attempt to update pipeline info including null values should fail
    assertThrows(
        ConstraintViolationException.class,
        () ->
            pipelinesService.adminUpdatePipelineWorkspace(
                pipelinesEnum, version, null, null, newWorkspaceName, null));
  }
}
