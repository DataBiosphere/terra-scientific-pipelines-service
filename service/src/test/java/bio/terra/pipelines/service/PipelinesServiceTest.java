package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.params.provider.Arguments.arguments;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineRuntimeMetadata;
import bio.terra.pipelines.db.repositories.PipelineRuntimeMetadataRepository;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.model.Pipeline;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
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
  @Autowired PipelineRuntimeMetadataRepository pipelineRuntimeMetadataRepository;
  @MockitoBean SamService samService;
  @MockitoBean RawlsService rawlsService;
  Integer arrayImputationNonHiddenInLiquiBasePipelineVersion = 1;
  Integer savedPipelineVersion = 100;

  @Test
  void getCorrectNumberOfPipelines() {
    // YAML defines one visible pipeline (v1) and one hidden pipeline (v2) by default.
    List<Pipeline> pipelineList = pipelinesService.getPipelines(false);
    assertEquals(1, pipelineList.size());
    String workspaceBillingProject = "testTerraProject";
    String workspaceName = "testTerraWorkspaceName";
    String workspaceStorageContainerName = "testWorkspaceStorageContainerUrl";
    String workspaceGoogleProject = "testWorkspaceGoogleProject";

    // save a DB-only version that is not present in YAML; it should not affect available pipelines
    pipelineRuntimeMetadataRepository.save(
        new PipelineRuntimeMetadata(
            PipelinesEnum.ARRAY_IMPUTATION,
            savedPipelineVersion,
            false,
            "1.2.1",
            workspaceBillingProject,
            workspaceName,
            workspaceStorageContainerName,
            workspaceGoogleProject));

    pipelineList = pipelinesService.getPipelines(false);
    assertEquals(1, pipelineList.size());
    Pipeline visiblePipeline = pipelineList.get(0);
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, visiblePipeline.getName());
    assertEquals(arrayImputationNonHiddenInLiquiBasePipelineVersion, visiblePipeline.getVersion());

    // test how many hidden pipelines exist
    pipelineList = pipelinesService.getPipelines(true);
    assertEquals(2, pipelineList.size());

    // save another DB-only hidden version; it should also not affect the YAML-defined list
    pipelineRuntimeMetadataRepository.save(
        new PipelineRuntimeMetadata(
            PipelinesEnum.ARRAY_IMPUTATION,
            savedPipelineVersion + 1,
            true,
            "1.2.1",
            workspaceBillingProject,
            workspaceName,
            workspaceStorageContainerName,
            workspaceGoogleProject));

    // make sure hidden pipeline is not returned when showHidden is false
    pipelineList = pipelinesService.getPipelines(false);
    assertEquals(1, pipelineList.size());

    pipelineList = pipelinesService.getPipelines(true);
    assertEquals(2, pipelineList.size());
  }

  @Test
  void getPipelineByKey() {
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;
    Pipeline p = pipelinesService.getPipeline(imputationPipeline, null, false);
    String pipelineKey = p.getKey();

    Pipeline pByKey = pipelinesService.getPipelineByKey(pipelineKey);
    assertEquals(p, pByKey);

    assertThrows(
        NotFoundException.class,
        () -> pipelinesService.getPipelineByKey("not_a_valid_pipeline_id"));
  }

  @Test
  void getPipelineUsesYamlMetadataAndDbMutableFields() {
    PipelineRuntimeMetadata dbPipeline =
        pipelineRuntimeMetadataRepository.findByNameAndVersion(PipelinesEnum.ARRAY_IMPUTATION, 1);
    dbPipeline.setToolVersion("9.9.9");
    dbPipeline.setWorkspaceBillingProject("db-billing-project");
    dbPipeline.setWorkspaceName("db-workspace-name");
    dbPipeline.setWorkspaceStorageContainerName("db-storage-container");
    dbPipeline.setWorkspaceGoogleProject("db-google-project");
    pipelineRuntimeMetadataRepository.save(dbPipeline);

    Pipeline hydratedPipeline =
        pipelinesService.getPipeline(PipelinesEnum.ARRAY_IMPUTATION, 1, true);

    assertEquals("Test Array Imputation", hydratedPipeline.getDisplayName());
    assertEquals("Test description", hydratedPipeline.getDescription());
    assertEquals("imputation", hydratedPipeline.getPipelineType());
    assertEquals("ImputationBeagle", hydratedPipeline.getToolName());

    assertEquals("9.9.9", hydratedPipeline.getToolVersion());
    assertEquals("db-billing-project", hydratedPipeline.getWorkspaceBillingProject());
    assertEquals("db-workspace-name", hydratedPipeline.getWorkspaceName());
    assertEquals("db-storage-container", hydratedPipeline.getWorkspaceStorageContainerName());
    assertEquals("db-google-project", hydratedPipeline.getWorkspaceGoogleProject());

    assertEquals(1, hydratedPipeline.getPipelineInputDefinitions().size());
    assertEquals(1, hydratedPipeline.getPipelineOutputDefinitions().size());
    assertTrue(
        hydratedPipeline.getPipelineInputDefinitions().stream()
            .anyMatch(input -> input.getName().equals("multiSampleVcf")));
    assertTrue(
        hydratedPipeline.getPipelineOutputDefinitions().stream()
            .anyMatch(output -> output.getName().equals("imputedMultiSampleVcf")));
  }

  @Test
  void getPipelineExistsFromYamlWithoutDbRow() {
    PipelineRuntimeMetadata versionTwoDbPipeline =
        pipelineRuntimeMetadataRepository.findByNameAndVersion(PipelinesEnum.ARRAY_IMPUTATION, 2);
    pipelineRuntimeMetadataRepository.deleteById(versionTwoDbPipeline.getId());

    Pipeline pipeline = pipelinesService.getPipeline(PipelinesEnum.ARRAY_IMPUTATION, 2, true);

    assertNull(pipeline.getKey());
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, pipeline.getName());
    assertEquals(2, pipeline.getVersion());
    assertEquals("All of Us + AnVIL Array Imputation", pipeline.getDisplayName());
    assertNull(pipeline.getToolVersion());
    assertNull(pipeline.getWorkspaceBillingProject());
    assertEquals(8, pipeline.getPipelineInputDefinitions().size());
    assertEquals(6, pipeline.getPipelineOutputDefinitions().size());
  }

  @Test
  void getHiddenPipeline() {
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;

    // version 2 is hidden in DB and should not be returned to non-admin callers
    assertThrows(
        NotFoundException.class, () -> pipelinesService.getPipeline(imputationPipeline, 2, false));

    Pipeline p = pipelinesService.getPipeline(imputationPipeline, 2, true);
    assertEquals(2, p.getVersion());

    // should grab non-hidden pipeline when pipeline version is not provided even if showHidden is
    // true
    p = pipelinesService.getPipeline(imputationPipeline, null, true);
    assertEquals(arrayImputationNonHiddenInLiquiBasePipelineVersion, p.getVersion());
  }

  @Test
  void getPipelineByNameAndVersion() {
    // save a DB-only version of the same pipeline; YAML should remain the source of truth
    pipelineRuntimeMetadataRepository.save(
        new PipelineRuntimeMetadata(
            PipelinesEnum.ARRAY_IMPUTATION,
            savedPipelineVersion,
            false,
            "1.2.1",
            null,
            null,
            null,
            null));
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;
    // this should return the highest version of the pipeline
    Pipeline nullVersionPipeline = pipelinesService.getPipeline(imputationPipeline, null, false);
    assertEquals(
        arrayImputationNonHiddenInLiquiBasePipelineVersion, nullVersionPipeline.getVersion());

    // this should return the specific version of the pipeline that exists
    Pipeline specificVersionPipeline =
        pipelinesService.getPipeline(
            imputationPipeline, arrayImputationNonHiddenInLiquiBasePipelineVersion, false);
    assertEquals(
        arrayImputationNonHiddenInLiquiBasePipelineVersion, specificVersionPipeline.getVersion());

    // if asking for unknown version pipeline combo, throw exception
    assertThrows(
        NotFoundException.class,
        () -> pipelinesService.getPipeline(imputationPipeline, savedPipelineVersion, false));
  }

  @Test
  void getLatestPipeline() {
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;
    // this should return the highest version of the pipeline
    Pipeline getLatestPipeline = pipelinesService.getLatestPipeline(imputationPipeline);
    assertEquals(
        arrayImputationNonHiddenInLiquiBasePipelineVersion, getLatestPipeline.getVersion());

    // a DB-only version should not affect the latest configured pipeline version
    pipelineRuntimeMetadataRepository.save(
        new PipelineRuntimeMetadata(
            PipelinesEnum.ARRAY_IMPUTATION, 100, false, "1.2.1", null, null, null, null));
    getLatestPipeline = pipelinesService.getLatestPipeline(imputationPipeline);
    assertEquals(1, getLatestPipeline.getVersion());

    PipelineRuntimeMetadata versionOne =
        pipelineRuntimeMetadataRepository.findByNameAndVersion(imputationPipeline, 1);
    PipelineRuntimeMetadata versionTwo =
        pipelineRuntimeMetadataRepository.findByNameAndVersion(imputationPipeline, 2);
    versionOne.setHidden(true);
    versionTwo.setHidden(false);
    pipelineRuntimeMetadataRepository.save(versionOne);
    pipelineRuntimeMetadataRepository.save(versionTwo);

    getLatestPipeline = pipelinesService.getLatestPipeline(imputationPipeline);
    assertEquals(2, getLatestPipeline.getVersion());
  }

  @Test
  void testPipelineToString() {
    // test .ToString() method on Pipeline Entity
    List<Pipeline> pipelineList = pipelinesService.getPipelines(false);
    assertEquals(1, pipelineList.size());
    for (Pipeline p : pipelineList) {
      assertEquals(
          String.format(
              "Pipeline[pipelineName=%s, version=%s, hidden=%s, displayName=%s, description=%s, pipelineType=%s, toolName=%s, toolVersion=%s, workspaceBillingProject=%s, workspaceName=%s, workspaceStorageContainerName=%s, workspaceGoogleProject=%s]",
              p.getName(),
              p.getVersion(),
              p.isHidden(),
              p.getDisplayName(),
              p.getDescription(),
              p.getPipelineType(),
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
              .append(p.getKey())
              .append(p.getName())
              .append(p.getVersion())
              .append(p.isHidden())
              .append(p.getDisplayName())
              .append(p.getDescription())
              .append(p.getPipelineType())
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

  @Test
  void getPipelinesOrderedByNameAndVersion() {
    pipelineRuntimeMetadataRepository.save(
        TestUtils.createTestPipelineRuntime(PipelinesEnum.ARRAY_IMPUTATION, 5, false, "1.2.4"));

    List<Pipeline> pipelineList = pipelinesService.getPipelines(false);
    assertEquals(1, pipelineList.size());

    List<Integer> actualVersions = pipelineList.stream().map(Pipeline::getVersion).toList();
    assertEquals(List.of(1), actualVersions);
  }

  @Test
  void getPipelinesOrderedByNameAndVersionIncludingHidden() {
    pipelineRuntimeMetadataRepository.save(
        TestUtils.createTestPipelineRuntime(PipelinesEnum.ARRAY_IMPUTATION, 5, true, "1.3.0"));

    List<Pipeline> pipelineList = pipelinesService.getPipelines(true);

    assertEquals(2, pipelineList.size());

    List<Integer> actualVersions = pipelineList.stream().map(Pipeline::getVersion).toList();
    assertEquals(List.of(2, 1), actualVersions);
  }
}
