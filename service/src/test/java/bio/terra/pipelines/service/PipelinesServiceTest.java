package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.params.provider.Arguments.arguments;
import static org.mockito.Mockito.when;

import bio.terra.cbas.model.OutputDestination;
import bio.terra.cbas.model.ParameterDefinition;
import bio.terra.cbas.model.ParameterTypeDefinition;
import bio.terra.cbas.model.ParameterTypeDefinitionArray;
import bio.terra.cbas.model.ParameterTypeDefinitionPrimitive;
import bio.terra.cbas.model.PrimitiveParameterValueType;
import bio.terra.cbas.model.WorkflowInputDefinition;
import bio.terra.cbas.model.WorkflowOutputDefinition;
import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.db.repositories.PipelinesRepository;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.rawls.model.WorkspaceDetails;
import jakarta.validation.ConstraintViolationException;
import java.util.ArrayList;
import java.util.List;
import java.util.UUID;
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

  @Test
  void getCorrectNumberOfPipelines() {
    // migrations insert one pipeline (imputation) so make sure we find it
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    assertEquals(1, pipelineList.size());
    UUID workspaceId = UUID.randomUUID();
    String workspaceBillingProject = "testTerraProject";
    String workspaceName = "testTerraWorkspaceName";
    String workspaceStorageContainerName = "testWorkspaceStorageContainerUrl";
    String workspaceGoogleProject = "testWorkspaceGoogleProject";

    // save a new version of the same pipeline
    pipelinesRepository.save(
        new Pipeline(
            PipelinesEnum.ARRAY_IMPUTATION,
            1,
            "pipelineDisplayName",
            "description",
            "pipelineType",
            "wdlUrl",
            "toolName",
            "1.2.1",
            workspaceId,
            workspaceBillingProject,
            workspaceName,
            workspaceStorageContainerName,
            workspaceGoogleProject,
            null,
            null));

    pipelineList = pipelinesService.getPipelines();
    assertEquals(2, pipelineList.size());
    Pipeline savedPipeline = pipelineList.get(1);
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, savedPipeline.getName());
    assertEquals(1, savedPipeline.getVersion());
    assertEquals("pipelineDisplayName", savedPipeline.getDisplayName());
    assertEquals("description", savedPipeline.getDescription());
    assertEquals("pipelineType", savedPipeline.getPipelineType());
    assertEquals("wdlUrl", savedPipeline.getWdlUrl());
    assertEquals("toolName", savedPipeline.getToolName());
    assertEquals("1.2.1", savedPipeline.getToolVersion());
    assertEquals(workspaceId, savedPipeline.getWorkspaceId());
    assertEquals(workspaceBillingProject, savedPipeline.getWorkspaceBillingProject());
    assertEquals(workspaceName, savedPipeline.getWorkspaceName());
    assertEquals(workspaceStorageContainerName, savedPipeline.getWorkspaceStorageContainerName());
    assertEquals(workspaceGoogleProject, savedPipeline.getWorkspaceGoogleProject());
  }

  @Test
  void getPipelineById() {
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;
    Pipeline p = pipelinesService.getPipeline(imputationPipeline, null);
    Long pipelineId = p.getId();

    Pipeline pById = pipelinesService.getPipelineById(pipelineId);
    assertEquals(p, pById);

    assertThrows(NotFoundException.class, () -> pipelinesService.getPipelineById(999L));
  }

  @Test
  void getPipelineByNameAndVersion() {
    // save a new version of the same pipeline that exists in the table
    pipelinesRepository.save(
        new Pipeline(
            PipelinesEnum.ARRAY_IMPUTATION,
            1,
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
            null,
            null));
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;
    // this should return the highest version of the pipeline
    Pipeline nullVersionPipeline = pipelinesService.getPipeline(imputationPipeline, null);
    assertEquals(1, nullVersionPipeline.getVersion());

    // this should return the specific version of the pipeline that exists
    Pipeline specificVersionPipeline = pipelinesService.getPipeline(imputationPipeline, 0);
    assertEquals(0, specificVersionPipeline.getVersion());

    // if asking for unknown version pipeline combo, throw exception
    assertThrows(
        NotFoundException.class, () -> pipelinesService.getPipeline(imputationPipeline, 999));
  }

  @Test
  void getLatestPipeline() {
    PipelinesEnum imputationPipeline = PipelinesEnum.ARRAY_IMPUTATION;
    // this should return the highest version of the pipeline
    Pipeline getLatestPipeline = pipelinesService.getLatestPipeline(imputationPipeline);
    assertEquals(0, getLatestPipeline.getVersion());

    // save a new version of the same pipeline that exists in the table
    pipelinesRepository.save(
        new Pipeline(
            PipelinesEnum.ARRAY_IMPUTATION,
            100,
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
            null,
            null));

    // this should return the new highest version of the pipeline
    getLatestPipeline = pipelinesService.getLatestPipeline(imputationPipeline);
    assertEquals(100, getLatestPipeline.getVersion());
  }

  @Test
  void testPipelineToString() {
    // test .ToString() method on Pipeline Entity
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    assertEquals(1, pipelineList.size());
    for (Pipeline p : pipelineList) {
      assertEquals(
          String.format(
              "Pipeline[pipelineName=%s, version=%s, displayName=%s, description=%s, pipelineType=%s, wdlUrl=%s, toolName=%s, toolVersion=%s, workspaceId=%s, workspaceBillingProject=%s, workspaceName=%s, workspaceStorageContainerName=%s, workspaceGoogleProject=%s]",
              p.getName(),
              p.getVersion(),
              p.getDisplayName(),
              p.getDescription(),
              p.getPipelineType(),
              p.getWdlUrl(),
              p.getToolName(),
              p.getToolVersion(),
              p.getWorkspaceId(),
              p.getWorkspaceBillingProject(),
              p.getWorkspaceName(),
              p.getWorkspaceStorageContainerName(),
              p.getWorkspaceGoogleProject()),
          p.toString());
    }
  }

  @Test
  void testPipelineHashCode() {
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    assertEquals(1, pipelineList.size());
    for (Pipeline p : pipelineList) {
      // 17 and 31 are hardcoded in this hashCode method of this class
      assertEquals(
          new HashCodeBuilder(17, 31)
              .append(p.getId())
              .append(p.getName())
              .append(p.getVersion())
              .append(p.getDisplayName())
              .append(p.getDescription())
              .append(p.getPipelineType())
              .append(p.getWdlUrl())
              .append(p.getToolName())
              .append(p.getToolVersion())
              .append(p.getWorkspaceId())
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
    Pipeline p = pipelinesService.getPipeline(pipelinesEnum, null);
    String newWorkspaceBillingProject = "newTestTerraProject";
    String newWorkspaceName = "newTestTerraWorkspaceName";
    String newToolVersion = "0.13.1";

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

    // mocks
    when(samService.getTeaspoonsServiceAccountToken()).thenReturn("fakeToken");
    when(rawlsService.getWorkspaceDetails(
            "fakeToken", newWorkspaceBillingProject, newWorkspaceName))
        .thenReturn(workspaceDetails);
    when(rawlsService.getWorkspaceBucketName(workspaceDetails))
        .thenReturn(newWorkspaceStorageContainerName);
    when(rawlsService.getWorkspaceGoogleProject(workspaceDetails))
        .thenReturn(newWorkspaceGoogleProject);

    // update pipeline workspace id
    pipelinesService.adminUpdatePipelineWorkspace(
        pipelinesEnum,
        p.getVersion(),
        newWorkspaceBillingProject,
        newWorkspaceName,
        newToolVersion);
    p = pipelinesService.getPipeline(pipelinesEnum, p.getVersion());

    // assert the workspace info has been updated
    assertEquals(newWorkspaceBillingProject, p.getWorkspaceBillingProject());
    assertEquals(newWorkspaceName, p.getWorkspaceName());
    assertEquals(newToolVersion, p.getToolVersion());

    // assert the fetched info from Rawls has been written
    assertEquals(newWorkspaceStorageContainerName, p.getWorkspaceStorageContainerName());
    assertEquals(newWorkspaceGoogleProject, p.getWorkspaceGoogleProject());
  }

  private static Stream<Arguments> badToolVersions() {
    return Stream.of(
        arguments("1.13.1"), // current imputation beagle pipeline is version 0 so this should fail
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
    Pipeline p = pipelinesService.getPipeline(pipelinesEnum, null);
    String newWorkspaceBillingProject = "newTestTerraProject";
    String newWorkspaceName = "newTestTerraWorkspaceName";
    int version = p.getVersion();
    assertThrows(
        ValidationException.class,
        () ->
            pipelinesService.adminUpdatePipelineWorkspace(
                pipelinesEnum,
                version,
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
        arguments("stringv0.1.32"));
  }

  @ParameterizedTest
  @MethodSource("goodToolVersions")
  void adminUpdatePipelineWorkspaceGoodToolVersion(String goodToolVersion) {
    PipelinesEnum pipelinesEnum = PipelinesEnum.ARRAY_IMPUTATION;
    Pipeline p = pipelinesService.getPipeline(pipelinesEnum, null);
    String newWorkspaceBillingProject = "newTestTerraProject";
    String newWorkspaceName = "newTestTerraWorkspaceName";
    assertDoesNotThrow(
        () ->
            pipelinesService.adminUpdatePipelineWorkspace(
                pipelinesEnum,
                p.getVersion(),
                newWorkspaceBillingProject,
                newWorkspaceName,
                goodToolVersion));
  }

  @Test
  void adminUpdatePipelineWorkspaceNullValueThrows() {
    PipelinesEnum pipelinesEnum = PipelinesEnum.ARRAY_IMPUTATION;
    Pipeline p = pipelinesService.getPipeline(pipelinesEnum, null);

    String newWorkspaceName = "newTestTerraWorkspaceName";
    int version = p.getVersion();

    // make sure the current pipeline does not have the workspace info we're trying to update with
    assertNotEquals(newWorkspaceName, p.getWorkspaceName());

    // attempt to update pipeline info including null values should fail
    assertThrows(
        ConstraintViolationException.class,
        () ->
            pipelinesService.adminUpdatePipelineWorkspace(
                pipelinesEnum, version, null, newWorkspaceName, null));
  }

  private static Stream<Arguments> mapVariableTypeToCbasParameterTypeArguments() {
    ParameterTypeDefinition stringParameterTypeResponse =
        new ParameterTypeDefinitionPrimitive()
            .primitiveType(PrimitiveParameterValueType.STRING)
            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE);
    ParameterTypeDefinition fileParameterTypeResponse =
        new ParameterTypeDefinitionPrimitive()
            .primitiveType(PrimitiveParameterValueType.FILE)
            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE);
    ParameterTypeDefinition integerParameterTypeResponse =
        new ParameterTypeDefinitionPrimitive()
            .primitiveType(PrimitiveParameterValueType.INT)
            .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE);
    ParameterTypeDefinition stringArrayParameterTypeResponse =
        new ParameterTypeDefinitionArray()
            .nonEmpty(true)
            .arrayType(
                new ParameterTypeDefinitionPrimitive()
                    .primitiveType(PrimitiveParameterValueType.STRING)
                    .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
            .type(ParameterTypeDefinition.TypeEnum.ARRAY);
    ParameterTypeDefinition fileArrayParameterTypeResponse =
        new ParameterTypeDefinitionArray()
            .nonEmpty(true)
            .arrayType(
                new ParameterTypeDefinitionPrimitive()
                    .primitiveType(PrimitiveParameterValueType.FILE)
                    .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
            .type(ParameterTypeDefinition.TypeEnum.ARRAY);
    return Stream.of(
        // arguments: type specification, expected response
        arguments(PipelineVariableTypesEnum.STRING, stringParameterTypeResponse),
        arguments(PipelineVariableTypesEnum.FILE, fileParameterTypeResponse),
        arguments(PipelineVariableTypesEnum.INTEGER, integerParameterTypeResponse),
        arguments(PipelineVariableTypesEnum.STRING_ARRAY, stringArrayParameterTypeResponse),
        arguments(PipelineVariableTypesEnum.FILE_ARRAY, fileArrayParameterTypeResponse));
  }

  @ParameterizedTest
  @MethodSource("mapVariableTypeToCbasParameterTypeArguments")
  void mapVariableTypeToCbasParameterType(
      PipelineVariableTypesEnum inputType, ParameterTypeDefinition expectedResponse) {
    assertEquals(expectedResponse, pipelinesService.mapVariableTypeToCbasParameterType(inputType));
  }

  @Test
  void prepareCbasWorkflowInputRecordLookupDefinitions() {
    List<PipelineInputDefinition> inputDefinitions = new ArrayList<>();
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L,
            "input1",
            "input_1",
            PipelineVariableTypesEnum.STRING,
            null,
            true,
            true,
            false,
            null));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L,
            "input2",
            "input_2",
            PipelineVariableTypesEnum.INTEGER,
            null,
            false,
            true,
            false,
            "1"));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L,
            "input3",
            "input_3",
            PipelineVariableTypesEnum.STRING_ARRAY,
            null,
            true,
            false,
            false,
            "[\"1\", \"2\"]"));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L,
            "input4",
            "input_4",
            PipelineVariableTypesEnum.FILE,
            ".vcf.gz",
            false,
            false,
            false,
            "fake/file.vcf.gz"));
    inputDefinitions.add(
        new PipelineInputDefinition(
            1L,
            "input5",
            "input_5",
            PipelineVariableTypesEnum.FILE_ARRAY,
            ".vcf.gz",
            true,
            true,
            false,
            null));

    String testWdlName = "aFakeWdl";

    List<WorkflowInputDefinition> cbasWorkflowInputDefinitions =
        pipelinesService.prepareCbasWorkflowInputRecordLookupDefinitions(
            inputDefinitions, testWdlName);

    assertEquals(inputDefinitions.size(), cbasWorkflowInputDefinitions.size());
    for (int i = 0; i < inputDefinitions.size(); i++) {
      // the input type should be the object returned by the mapVariableTypeToCbasParameterType
      // method
      assertEquals(
          pipelinesService.mapVariableTypeToCbasParameterType(inputDefinitions.get(i).getType()),
          cbasWorkflowInputDefinitions.get(i).getInputType());
      // the input name should be the wdl name concatenated with the input name
      assertEquals(
          "%s.%s".formatted(testWdlName, inputDefinitions.get(i).getWdlVariableName()),
          cbasWorkflowInputDefinitions.get(i).getInputName());
      // the source should be a record lookup
      assertEquals(
          ParameterDefinition.TypeEnum.RECORD_LOOKUP,
          cbasWorkflowInputDefinitions.get(i).getSource().getType());
    }
  }

  @Test
  void prepareCbasWorkflowOutputRecordUpdateDefinitions() {
    List<PipelineOutputDefinition> outputDefinitions = new ArrayList<>();
    outputDefinitions.add(
        new PipelineOutputDefinition(1L, "output1", "output_1", PipelineVariableTypesEnum.FILE));
    outputDefinitions.add(
        new PipelineOutputDefinition(1L, "output2", "output_2", PipelineVariableTypesEnum.STRING));
    String testWdlName = "aFakeWdl";

    List<WorkflowOutputDefinition> cbasWorkflowOutputDefinitions =
        pipelinesService.prepareCbasWorkflowOutputRecordUpdateDefinitions(
            outputDefinitions, testWdlName);

    assertEquals(outputDefinitions.size(), cbasWorkflowOutputDefinitions.size());
    for (int i = 0; i < outputDefinitions.size(); i++) {
      // the output name should be the wdl name concatenated with the output name
      assertEquals(
          "%s.%s".formatted(testWdlName, outputDefinitions.get(i).getWdlVariableName()),
          cbasWorkflowOutputDefinitions.get(i).getOutputName());
      // the output type should be the object returned by the mapVariableTypeToCbasParameterType
      // method
      assertEquals(
          pipelinesService.mapVariableTypeToCbasParameterType(outputDefinitions.get(i).getType()),
          cbasWorkflowOutputDefinitions.get(i).getOutputType());
      // the destination type should be a record update
      assertEquals(
          OutputDestination.TypeEnum.RECORD_UPDATE,
          cbasWorkflowOutputDefinitions.get(i).getDestination().getType());
    }
  }
}
