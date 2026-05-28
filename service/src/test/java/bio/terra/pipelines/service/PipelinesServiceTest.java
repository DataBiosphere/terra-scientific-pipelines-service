package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.params.provider.Arguments.arguments;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineRuntimeMetadata;
import bio.terra.pipelines.db.repositories.PipelineRuntimeMetadataRepository;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.rawls.model.WorkspaceDetails;
import jakarta.validation.ConstraintViolationException;
import java.util.List;
import java.util.stream.Stream;
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

  // Pipeline YAML configuration in tests is loaded only from test/resources/pipelines-config.yml:
  //   array_imputation v0 and v1; low_pass_imputation v1 → 3 total YAML versions.
  //
  // After Liquibase migrations, pipeline_runtime_metadata is seeded with:
  //   array_imputation_v1 (hidden=false), array_imputation_v2 (hidden=true),
  //   low_pass_imputation_v1 (hidden=true)
  //   array_imputation_v0 has NO runtime metadata row → defaults to hidden=false
  // Runtime metadata rows without a matching YAML definition (e.g. array_imputation_v2) are
  // ignored.
  //
  // Initial visible (showHidden=false): ai_v0, ai_v1 (2 non-hidden)
  // Initial hidden: lpi_v1 (hidden=true)
  // Total (showHidden=true): 3
  static final int ARRAY_IMP_VERSION_0 = 0;
  static final int ARRAY_IMP_VERSION_1 = 1;
  static final int LOW_PASS_IMP_VERSION_1 = 1;
  static final String ARRAY_IMP_V0_KEY = "array_imputation_v0";
  static final String ARRAY_IMP_V1_KEY = "array_imputation_v1";
  static final String LOW_PASS_IMP_V1_KEY = "low_pass_imputation_v1";

  /**
   * Upsert a PipelineRuntimeMetadata row for {@code pipelineKey} with the given hidden flag. Other
   * runtime fields (workspace details, toolVersion) are left null unless the caller sets them on
   * the returned object separately.
   */
  private void saveRuntimeMetadata(String pipelineKey, boolean hidden) {
    PipelineRuntimeMetadata meta = new PipelineRuntimeMetadata();
    meta.setPipelineKey(pipelineKey);
    meta.setHidden(hidden);
    pipelineRuntimeMetadataRepository.save(meta);
  }

  @Test
  void getCorrectNumberOfPipelines() {
    // Initial state: 3 YAML-defined versions total; 2 are non-hidden (ai_v0, ai_v1)
    List<Pipeline> pipelineListIncludingHidden = pipelinesService.getPipelines(true);
    assertEquals(3, pipelineListIncludingHidden.size());

    List<Pipeline> pipelineList = pipelinesService.getPipelines(false);
    assertEquals(2, pipelineList.size());

    // verify ordering: higher version first within the same pipeline name
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, pipelineList.get(0).getName());
    assertEquals(ARRAY_IMP_VERSION_1, pipelineList.get(0).getVersion());
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, pipelineList.get(1).getName());
    assertEquals(ARRAY_IMP_VERSION_0, pipelineList.get(1).getVersion());

    // save runtime metadata with workspace details for ai_v0 and verify they are visible
    String workspaceBillingProject = "testTerraProject";
    String workspaceName = "testTerraWorkspaceName";
    String workspaceStorageContainerName = "testWorkspaceStorageContainerUrl";
    String workspaceGoogleProject = "testWorkspaceGoogleProject";
    String toolVersion = "1.2.1";
    PipelineRuntimeMetadata meta = new PipelineRuntimeMetadata();
    meta.setPipelineKey(ARRAY_IMP_V0_KEY);
    meta.setHidden(false);
    meta.setToolVersion(toolVersion);
    meta.setWorkspaceBillingProject(workspaceBillingProject);
    meta.setWorkspaceName(workspaceName);
    meta.setWorkspaceStorageContainerName(workspaceStorageContainerName);
    meta.setWorkspaceGoogleProject(workspaceGoogleProject);
    pipelineRuntimeMetadataRepository.save(meta);

    Pipeline ai_v0 =
        pipelinesService.getPipeline(PipelinesEnum.ARRAY_IMPUTATION, ARRAY_IMP_VERSION_0, false);
    assertEquals(workspaceBillingProject, ai_v0.getWorkspaceBillingProject());
    assertEquals(workspaceName, ai_v0.getWorkspaceName());
    assertEquals(workspaceStorageContainerName, ai_v0.getWorkspaceStorageContainerName());
    assertEquals(workspaceGoogleProject, ai_v0.getWorkspaceGoogleProject());
    assertEquals(toolVersion, ai_v0.getToolVersion());

    // make lpi_v1 non-hidden → now 3 are visible (ai_v0, ai_v1, lpi_v1)
    saveRuntimeMetadata(LOW_PASS_IMP_V1_KEY, false);
    pipelineList = pipelinesService.getPipelines(false);
    assertEquals(3, pipelineList.size());

    // make ai_v0 hidden → back to 2 non-hidden (ai_v1 and lpi_v1)
    saveRuntimeMetadata(ARRAY_IMP_V0_KEY, true);
    pipelineList = pipelinesService.getPipelines(false);
    assertEquals(2, pipelineList.size());

    // total with showHidden=true is still 3
    pipelineList = pipelinesService.getPipelines(true);
    assertEquals(3, pipelineList.size());
  }

  @Test
  void getPipelineById() {
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;
    // getPipeline(null, false) returns ai_v1 (highest non-hidden version)
    Pipeline p = pipelinesService.getPipeline(imputationPipeline, null, false);
    Long pipelineId = p.getId();

    Pipeline pById = pipelinesService.getPipelineById(pipelineId);
    assertEquals(p, pById);

    assertThrows(NotFoundException.class, () -> pipelinesService.getPipelineById(999L));
  }

  @Test
  void getHiddenPipeline() {
    // ai_v1 is initially non-hidden; make it hidden via runtime metadata
    saveRuntimeMetadata(ARRAY_IMP_V1_KEY, true);

    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;

    // should not be able to grab ai_v1 when showHidden is false
    assertThrows(
        NotFoundException.class,
        () -> pipelinesService.getPipeline(imputationPipeline, ARRAY_IMP_VERSION_1, false));

    // should be able to grab ai_v1 with showHidden=true
    Pipeline p = pipelinesService.getPipeline(imputationPipeline, ARRAY_IMP_VERSION_1, true);
    assertEquals(ARRAY_IMP_VERSION_1, p.getVersion());

    // getPipeline(null, false) falls back to ai_v0 (still visible; no runtime metadata → not
    // hidden)
    p = pipelinesService.getPipeline(imputationPipeline, null, false);
    assertEquals(ARRAY_IMP_VERSION_0, p.getVersion());

    // getPipeline(null, true) returns the highest version regardless of hidden status (= v1)
    p = pipelinesService.getPipeline(imputationPipeline, null, true);
    assertEquals(ARRAY_IMP_VERSION_1, p.getVersion());
  }

  @Test
  void getPipelineByNameAndVersion() {
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;

    // ai_v0 and ai_v1 are both non-hidden; highest non-hidden version is 1
    Pipeline nullVersionPipeline = pipelinesService.getPipeline(imputationPipeline, null, false);
    assertEquals(ARRAY_IMP_VERSION_1, nullVersionPipeline.getVersion());

    // specific version lookup for the lower version
    Pipeline specificVersionPipeline =
        pipelinesService.getPipeline(imputationPipeline, ARRAY_IMP_VERSION_0, false);
    assertEquals(ARRAY_IMP_VERSION_0, specificVersionPipeline.getVersion());

    // unknown version is rejected when no YAML pipeline definition exists
    assertThrows(
        IllegalArgumentException.class,
        () -> pipelinesService.getPipeline(imputationPipeline, 999, false));

    // make ai_v0 hidden; getPipeline(null, false) should still return ai_v1
    saveRuntimeMetadata(ARRAY_IMP_V0_KEY, true);
    nullVersionPipeline = pipelinesService.getPipeline(imputationPipeline, null, false);
    assertEquals(ARRAY_IMP_VERSION_1, nullVersionPipeline.getVersion());

    // getPipeline(v0, false) now throws because ai_v0 is hidden
    assertThrows(
        NotFoundException.class,
        () -> pipelinesService.getPipeline(imputationPipeline, ARRAY_IMP_VERSION_0, false));

    // getPipeline(v0, true) returns ai_v0 even though it is hidden
    Pipeline hiddenPipeline =
        pipelinesService.getPipeline(imputationPipeline, ARRAY_IMP_VERSION_0, true);
    assertEquals(ARRAY_IMP_VERSION_0, hiddenPipeline.getVersion());
  }

  @Test
  void getLatestPipeline() {
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;

    // getLatestPipeline uses showHidden=true; with test YAML having v0 and v1, returns v1
    Pipeline latestPipeline = pipelinesService.getLatestPipeline(imputationPipeline);
    assertEquals(ARRAY_IMP_VERSION_1, latestPipeline.getVersion());

    // make ai_v1 hidden; getLatestPipeline still returns v1 because showHidden=true
    saveRuntimeMetadata(ARRAY_IMP_V1_KEY, true);
    latestPipeline = pipelinesService.getLatestPipeline(imputationPipeline);
    assertEquals(ARRAY_IMP_VERSION_1, latestPipeline.getVersion());

    // getPipeline(null, false) now skips the hidden ai_v1 and returns ai_v0
    Pipeline latestVisible = pipelinesService.getPipeline(imputationPipeline, null, false);
    assertEquals(ARRAY_IMP_VERSION_0, latestVisible.getVersion());
  }

  @Test
  void adminUpdatePipelineWorkspace() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.ARRAY_IMPUTATION;
    // getPipeline(null, true) returns ai_v1 (highest version, showHidden=true includes hidden)
    Pipeline p = pipelinesService.getPipeline(pipelinesEnum, null, true);
    assertEquals(ARRAY_IMP_VERSION_1, p.getVersion());

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

    // set up mocks
    when(samService.getTeaspoonsServiceAccountToken()).thenReturn("fakeToken");
    when(rawlsService.getWorkspaceDetails(
            "fakeToken", newWorkspaceBillingProject, newWorkspaceName))
        .thenReturn(workspaceDetails);
    when(rawlsService.getWorkspaceBucketName(workspaceDetails))
        .thenReturn(newWorkspaceStorageContainerName);
    when(rawlsService.getWorkspaceGoogleProject(workspaceDetails))
        .thenReturn(newWorkspaceGoogleProject);

    // update pipeline values, including hiding ai_v1
    pipelinesService.adminUpdatePipelineWorkspace(
        pipelinesEnum,
        ARRAY_IMP_VERSION_1,
        newHidden,
        newWorkspaceBillingProject,
        newWorkspaceName,
        newToolVersion);
    p = pipelinesService.getPipeline(pipelinesEnum, ARRAY_IMP_VERSION_1, true);

    // assert workspace info has been updated
    assertEquals(newWorkspaceBillingProject, p.getWorkspaceBillingProject());
    assertEquals(newWorkspaceName, p.getWorkspaceName());
    assertEquals(newToolVersion, p.getToolVersion());
    assertEquals(newWorkspaceStorageContainerName, p.getWorkspaceStorageContainerName());
    assertEquals(newWorkspaceGoogleProject, p.getWorkspaceGoogleProject());
    assertTrue(p.isHidden());

    // ai_v1 is now hidden; getPipeline(null, false) should fall back to the visible ai_v0
    Pipeline visiblePipeline = pipelinesService.getPipeline(pipelinesEnum, null, false);
    assertEquals(ARRAY_IMP_VERSION_0, visiblePipeline.getVersion());

    // update with null isHidden value – should not change visibility
    pipelinesService.adminUpdatePipelineWorkspace(
        pipelinesEnum,
        ARRAY_IMP_VERSION_1,
        null,
        newWorkspaceBillingProject,
        newWorkspaceName,
        newToolVersion);
    p = pipelinesService.getPipeline(pipelinesEnum, ARRAY_IMP_VERSION_1, true);
    assertTrue(p.isHidden());

    // update with isHidden=false – ai_v1 becomes visible again
    pipelinesService.adminUpdatePipelineWorkspace(
        pipelinesEnum,
        ARRAY_IMP_VERSION_1,
        false,
        newWorkspaceBillingProject,
        newWorkspaceName,
        newToolVersion);
    p = pipelinesService.getPipeline(pipelinesEnum, ARRAY_IMP_VERSION_1, true);
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
        arguments("hiiv.1.4"),
        arguments(".."),
        arguments("..3"));
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
        arguments("1.13.1"),
        arguments("PipelineHasANumber8InIt_v1.2.3"));
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
    // Make lpi_v1 non-hidden so we can test cross-pipeline ordering
    saveRuntimeMetadata(LOW_PASS_IMP_V1_KEY, false);

    List<Pipeline> pipelineList = pipelinesService.getPipelines(false);
    assertEquals(3, pipelineList.size());

    // Verify array_imputation versions are in descending order
    List<Integer> arrayImputationVersions =
        pipelineList.stream()
            .filter(p -> p.getName().equals(PipelinesEnum.ARRAY_IMPUTATION))
            .map(Pipeline::getVersion)
            .toList();
    assertEquals(List.of(ARRAY_IMP_VERSION_1, ARRAY_IMP_VERSION_0), arrayImputationVersions);

    // Verify pipeline name ordering: ARRAY_IMPUTATION (alphabetically first) before
    // LOW_PASS_IMPUTATION
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, pipelineList.get(0).getName());
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, pipelineList.get(1).getName());
    assertEquals(PipelinesEnum.LOW_PASS_IMPUTATION, pipelineList.get(2).getName());
  }

  @Test
  void getPipelinesOrderedByNameAndVersionIncludingHidden() {
    // Test YAML defines 3 versions: ai_v0 and ai_v1 (non-hidden), lpi_v1 (hidden)
    List<Pipeline> pipelineList = pipelinesService.getPipelines(true);
    assertEquals(3, pipelineList.size());

    // Verify array_imputation versions are in descending order
    List<Integer> arrayImputationVersions =
        pipelineList.stream()
            .filter(p -> p.getName().equals(PipelinesEnum.ARRAY_IMPUTATION))
            .map(Pipeline::getVersion)
            .toList();
    assertEquals(List.of(ARRAY_IMP_VERSION_1, ARRAY_IMP_VERSION_0), arrayImputationVersions);

    // Verify low_pass_imputation appears after array_imputation (alphabetical ordering)
    List<Integer> lpiVersions =
        pipelineList.stream()
            .filter(p -> p.getName().equals(PipelinesEnum.LOW_PASS_IMPUTATION))
            .map(Pipeline::getVersion)
            .toList();
    assertEquals(List.of(LOW_PASS_IMP_VERSION_1), lpiVersions);
  }
}
