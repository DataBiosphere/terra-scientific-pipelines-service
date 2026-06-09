package bio.terra.pipelines.service;

import static bio.terra.pipelines.common.utils.PipelineKeyUtils.buildPipelineKey;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.model.Pipeline;
import bio.terra.pipelines.model.PipelineInputDefinition;
import bio.terra.pipelines.model.PipelineOutputDefinition;
import bio.terra.pipelines.stairway.steps.utils.ToolConfig;
import bio.terra.pipelines.testutils.BaseTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.math.BigDecimal;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

class ToolConfigServiceTest extends BaseTest {

  @Mock private PipelineConfigurations pipelineConfigurations;

  private ToolConfigService toolConfigService;

  private final PipelinesEnum pipelineName = PipelinesEnum.ARRAY_IMPUTATION;
  private final String toolName = TestUtils.TEST_TOOL_NAME_1;
  private final String toolVersion = TestUtils.TEST_TOOL_VERSION_1;
  private final List<PipelineInputDefinition> pipelineInputDefinitions =
      TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST;
  private final List<PipelineOutputDefinition> pipelineOutputDefinitions =
      TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST;
  private final String monitoringScriptPath = "gs://path/to/monitoring/script.sh";
  private final BigDecimal memoryRetryMultiplierPipeline = BigDecimal.valueOf(1.5);
  private final Long pollingIntervalSecondsMainTool = 30L;
  private final boolean useCallCachingMainTool = true;
  private final boolean deleteIntermediateFilesMainTool = false;
  private final Long pollingIntervalSecondsQuota = 5L;
  private final boolean useCallCachingQuota = false;
  private final Long pollingIntervalSecondsInputQc = 1L;
  private final boolean useCallCachingInputQc = true;

  private final int arrayImputationPipelineVersion = 2;
  private final int lowPassImputationPipelineVersion = 10;

  @BeforeEach
  void setUp() {
    MockitoAnnotations.openMocks(this);
    toolConfigService = new ToolConfigService(pipelineConfigurations);

    // mock array imputation config
    PipelineConfigurations.WdlBasedPipelineConfiguration arrayImputationPipelineConfiguration =
        buildPipelineConfigWithMemoryMultiplier(memoryRetryMultiplierPipeline);
    when(pipelineConfigurations.getPipelineConfiguration(
            "array_imputation_v%s".formatted(arrayImputationPipelineVersion)))
        .thenReturn(arrayImputationPipelineConfiguration);

    // mock low pass imputation config
    PipelineConfigurations.WdlBasedPipelineConfiguration lowPassImputationPipelineConfiguration =
        buildPipelineConfigWithMemoryMultiplier(memoryRetryMultiplierPipeline);
    when(pipelineConfigurations.getPipelineConfiguration(
            "low_pass_imputation_v%s".formatted(lowPassImputationPipelineVersion)))
        .thenReturn(lowPassImputationPipelineConfiguration);

    // mock pipelinesCommonConfiguration
    PipelineConfigurations.PipelinesCommonConfiguration pipelinesCommonConfiguration =
        mock(PipelineConfigurations.PipelinesCommonConfiguration.class);
    when(pipelineConfigurations.getCommon()).thenReturn(pipelinesCommonConfiguration);
    when(pipelinesCommonConfiguration.getQuotaConsumedPollingIntervalSeconds())
        .thenReturn(pollingIntervalSecondsQuota);
    when(pipelinesCommonConfiguration.isQuotaConsumedUseCallCaching())
        .thenReturn(useCallCachingQuota);
    when(pipelinesCommonConfiguration.getInputQcPollingIntervalSeconds())
        .thenReturn(pollingIntervalSecondsInputQc);
    when(pipelinesCommonConfiguration.isInputQcUseCallCaching()).thenReturn(useCallCachingInputQc);
    when(pipelinesCommonConfiguration.getMainToolPollingIntervalSeconds())
        .thenReturn(pollingIntervalSecondsMainTool);
    when(pipelinesCommonConfiguration.isMainToolUseCallCaching())
        .thenReturn(useCallCachingMainTool);
    when(pipelinesCommonConfiguration.isMainToolDeleteIntermediateFiles())
        .thenReturn(deleteIntermediateFilesMainTool);
    when(pipelinesCommonConfiguration.getMonitoringScriptPath()).thenReturn(monitoringScriptPath);
  }

  @Test
  void testGetPipelineMainToolConfigWithArrayImputationPipeline() {
    Pipeline pipeline =
        Pipeline.builder()
            .name(pipelineName)
            .version(arrayImputationPipelineVersion)
            .pipelineKey(buildPipelineKey(pipelineName, arrayImputationPipelineVersion))
            .toolName(toolName)
            .toolVersion(toolVersion)
            .inputDefinitions(pipelineInputDefinitions)
            .outputDefinitions(pipelineOutputDefinitions)
            .build();

    ToolConfig toolConfig = toolConfigService.getPipelineMainToolConfig(pipeline);

    assertEquals(toolName, toolConfig.methodName());
    assertEquals(toolVersion, toolConfig.methodVersion());
    assertEquals(
        "%s_v%s".formatted(toolName, arrayImputationPipelineVersion),
        toolConfig.methodNameWithPipelineVersion());
    assertEquals(
        "%s_v%s".formatted(pipelineName.getLowerCaseValue(), arrayImputationPipelineVersion),
        toolConfig.dataTableEntityName());
    assertEquals(pipelineInputDefinitions, toolConfig.inputDefinitions());
    assertEquals(pipelineOutputDefinitions, toolConfig.outputDefinitions());
    assertEquals(useCallCachingMainTool, toolConfig.callCache());
    assertEquals(monitoringScriptPath, toolConfig.monitoringScriptPath());
    assertEquals(deleteIntermediateFilesMainTool, toolConfig.deleteIntermediateOutputFiles());
    assertEquals(memoryRetryMultiplierPipeline, toolConfig.memoryRetryMultiplier());
    assertEquals(pollingIntervalSecondsMainTool, toolConfig.pollingIntervalSeconds());
  }

  @Test
  void testGetPipelineMainToolConfigWithLowPassImputationPipeline() {
    PipelinesEnum lowPassPipelineName = PipelinesEnum.LOW_PASS_IMPUTATION;

    Pipeline pipeline =
        Pipeline.builder()
            .name(lowPassPipelineName)
            .version(lowPassImputationPipelineVersion)
            .pipelineKey(buildPipelineKey(lowPassPipelineName, lowPassImputationPipelineVersion))
            .toolName(toolName)
            .toolVersion(toolVersion)
            .inputDefinitions(pipelineInputDefinitions)
            .outputDefinitions(pipelineOutputDefinitions)
            .build();

    ToolConfig toolConfig = toolConfigService.getPipelineMainToolConfig(pipeline);

    assertEquals(toolName, toolConfig.methodName());
    assertEquals(toolVersion, toolConfig.methodVersion());
    assertEquals(
        "%s_v%s".formatted(toolName, lowPassImputationPipelineVersion),
        toolConfig.methodNameWithPipelineVersion());
    assertEquals(
        "%s_v%s"
            .formatted(lowPassPipelineName.getLowerCaseValue(), lowPassImputationPipelineVersion),
        toolConfig.dataTableEntityName());
    assertEquals(pipelineInputDefinitions, toolConfig.inputDefinitions());
    assertEquals(pipelineOutputDefinitions, toolConfig.outputDefinitions());
    assertEquals(useCallCachingMainTool, toolConfig.callCache());
    assertEquals(monitoringScriptPath, toolConfig.monitoringScriptPath());
    assertEquals(deleteIntermediateFilesMainTool, toolConfig.deleteIntermediateOutputFiles());
    assertEquals(memoryRetryMultiplierPipeline, toolConfig.memoryRetryMultiplier());
    assertEquals(pollingIntervalSecondsMainTool, toolConfig.pollingIntervalSeconds());
  }

  @Test
  void testGetQuotaConsumedToolConfig() {
    // TEST_ARRAY_IMPUTATION_PIPELINE_1 has pipeline version 1
    Pipeline pipeline =
        TestUtils.TEST_ARRAY_IMPUTATION_PIPELINE_1.toBuilder().toolVersion(toolVersion).build();

    ToolConfig toolConfig = toolConfigService.getQuotaConsumedToolConfig(pipeline);

    List<PipelineOutputDefinition> expectedOutputDefinitions =
        List.of(
            PipelineOutputDefinition.builder()
                .name("quotaConsumed")
                .wdlVariableName("quota_consumed")
                .type(PipelineVariableTypesEnum.INTEGER)
                .isRequired(true)
                .build());
    assertEquals("QuotaConsumed", toolConfig.methodName());
    assertEquals(toolVersion, toolConfig.methodVersion());
    assertEquals(
        "QuotaConsumed_v%s".formatted(TestUtils.TEST_PIPELINE_VERSION_1),
        toolConfig.methodNameWithPipelineVersion());
    assertEquals(
        "array_imputation_v%s".formatted(TestUtils.TEST_PIPELINE_VERSION_1),
        toolConfig.dataTableEntityName());
    assertEquals(TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST, toolConfig.inputDefinitions());
    assertEquals(expectedOutputDefinitions, toolConfig.outputDefinitions());
    assertEquals(useCallCachingQuota, toolConfig.callCache());
    assertEquals(monitoringScriptPath, toolConfig.monitoringScriptPath());
    assertTrue(toolConfig.deleteIntermediateOutputFiles());
    assertNull(toolConfig.memoryRetryMultiplier());
    assertEquals(pollingIntervalSecondsQuota, toolConfig.pollingIntervalSeconds());
  }

  @Test
  void getInputQcToolConfig() {
    // TEST_ARRAY_IMPUTATION_PIPELINE_1 has pipeline version 1
    Pipeline pipeline =
        TestUtils.TEST_ARRAY_IMPUTATION_PIPELINE_1.toBuilder().toolVersion(toolVersion).build();

    ToolConfig toolConfig = toolConfigService.getInputQcToolConfig(pipeline);

    List<PipelineOutputDefinition> expectedOutputDefinitions =
        List.of(
            PipelineOutputDefinition.builder()
                .name("passesQc")
                .wdlVariableName("passes_qc")
                .type(PipelineVariableTypesEnum.BOOLEAN)
                .isRequired(true)
                .build(),
            PipelineOutputDefinition.builder()
                .name("qcMessages")
                .wdlVariableName("qc_messages")
                .type(PipelineVariableTypesEnum.STRING)
                .isRequired(false)
                .build());
    assertEquals("InputQC", toolConfig.methodName());
    assertEquals(toolVersion, toolConfig.methodVersion());
    assertEquals(
        "InputQC_v%s".formatted(TestUtils.TEST_PIPELINE_VERSION_1),
        toolConfig.methodNameWithPipelineVersion());
    assertEquals(
        "array_imputation_v%s".formatted(TestUtils.TEST_PIPELINE_VERSION_1),
        toolConfig.dataTableEntityName());
    assertEquals(TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST, toolConfig.inputDefinitions());
    assertEquals(expectedOutputDefinitions, toolConfig.outputDefinitions());
    assertEquals(useCallCachingInputQc, toolConfig.callCache());
    assertEquals(monitoringScriptPath, toolConfig.monitoringScriptPath());
    assertTrue(toolConfig.deleteIntermediateOutputFiles());
    assertNull(toolConfig.memoryRetryMultiplier());
    assertEquals(pollingIntervalSecondsInputQc, toolConfig.pollingIntervalSeconds());
  }

  private PipelineConfigurations.WdlBasedPipelineConfiguration
      buildPipelineConfigWithMemoryMultiplier(BigDecimal memoryRetryMultiplier) {
    PipelineConfigurations.WdlBasedPipelineConfiguration config =
        new PipelineConfigurations.WdlBasedPipelineConfiguration();
    config.setMemoryRetryMultiplier(memoryRetryMultiplier);
    return config;
  }
}
