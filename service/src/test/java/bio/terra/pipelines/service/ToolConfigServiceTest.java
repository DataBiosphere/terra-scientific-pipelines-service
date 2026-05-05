package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import bio.terra.pipelines.stairway.steps.utils.ToolConfig;
import bio.terra.pipelines.testutils.BaseTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.math.BigDecimal;
import java.util.List;
import java.util.Map;
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
  private final boolean useCallCachingPipeline = true;
  private final String monitoringScriptPath = "gs://path/to/monitoring/script.sh";
  private final boolean deleteIntermediateFilesPipeline = false;
  private final BigDecimal memoryRetryMultiplierPipeline = BigDecimal.valueOf(1.5);
  private final Long pollingIntervalSecondsPipeline = 30L;
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
    PipelineConfigurations.WdlBasedPipelineConfig arrayImputationPipelineConfig =
        new PipelineConfigurations.WdlBasedPipelineConfig(
            pollingIntervalSecondsPipeline,
            List.of(),
            "",
            Map.of(),
            useCallCachingPipeline,
            deleteIntermediateFilesPipeline,
            memoryRetryMultiplierPipeline);
    Map<String, PipelineConfigurations.WdlBasedPipelineConfig> arrayImputationConfigMap =
        Map.of(String.valueOf(arrayImputationPipelineVersion), arrayImputationPipelineConfig);
    when(pipelineConfigurations.getArrayImputation()).thenReturn(arrayImputationConfigMap);

    // mock low pass imputation config
    PipelineConfigurations.WdlBasedPipelineConfig lowPassImputationPipelineConfig =
        new PipelineConfigurations.WdlBasedPipelineConfig(
            pollingIntervalSecondsPipeline,
            List.of(),
            "",
            Map.of(),
            useCallCachingPipeline,
            deleteIntermediateFilesPipeline,
            memoryRetryMultiplierPipeline);
    Map<String, PipelineConfigurations.WdlBasedPipelineConfig> lowPassImputationConfigMap =
        Map.of(String.valueOf(lowPassImputationPipelineVersion), lowPassImputationPipelineConfig);
    when(pipelineConfigurations.getLowPassImputation()).thenReturn(lowPassImputationConfigMap);

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
    when(pipelinesCommonConfiguration.getMonitoringScriptPath()).thenReturn(monitoringScriptPath);
  }

  @Test
  void testGetPipelineMainToolConfigWithArrayImputationPipeline() {
    // create pipeline
    Pipeline pipeline = new Pipeline();
    pipeline.setName(pipelineName);
    pipeline.setVersion(arrayImputationPipelineVersion);
    pipeline.setToolName(toolName);
    pipeline.setToolVersion(toolVersion);
    pipeline.setPipelineInputDefinitions(pipelineInputDefinitions);
    pipeline.setPipelineOutputDefinitions(pipelineOutputDefinitions);

    // create main tool config
    ToolConfig toolConfig = toolConfigService.getPipelineMainToolConfig(pipeline);

    // check values
    assertEquals(toolName, toolConfig.methodName());
    assertEquals(toolVersion, toolConfig.methodVersion());
    assertEquals(
        "%s_v%s".formatted(toolName, arrayImputationPipelineVersion),
        toolConfig.methodNameWithPipelineVersion());
    assertEquals(
        "%s_v%s".formatted(pipelineName.getValue(), arrayImputationPipelineVersion),
        toolConfig.dataTableEntityName());
    assertEquals(pipelineInputDefinitions, toolConfig.inputDefinitions());
    assertEquals(pipelineOutputDefinitions, toolConfig.outputDefinitions());
    assertEquals(useCallCachingPipeline, toolConfig.callCache());
    assertEquals(monitoringScriptPath, toolConfig.monitoringScriptPath());
    assertEquals(deleteIntermediateFilesPipeline, toolConfig.deleteIntermediateOutputFiles());
    assertEquals(memoryRetryMultiplierPipeline, toolConfig.memoryRetryMultiplier());
    assertEquals(pollingIntervalSecondsPipeline, toolConfig.pollingIntervalSeconds());
  }

  @Test
  void testGetPipelineMainToolConfigWithLowPassImputationPipeline() {
    // create pipeline
    PipelinesEnum lowPassPipelineName = PipelinesEnum.LOW_PASS_IMPUTATION;

    Pipeline pipeline = new Pipeline();
    pipeline.setName(lowPassPipelineName);
    pipeline.setVersion(lowPassImputationPipelineVersion);
    pipeline.setToolName(toolName);
    pipeline.setToolVersion(toolVersion);
    pipeline.setPipelineInputDefinitions(pipelineInputDefinitions);
    pipeline.setPipelineOutputDefinitions(pipelineOutputDefinitions);

    // create main tool config
    ToolConfig toolConfig = toolConfigService.getPipelineMainToolConfig(pipeline);

    // check values
    assertEquals(toolName, toolConfig.methodName());
    assertEquals(toolVersion, toolConfig.methodVersion());
    assertEquals(
        "%s_v%s".formatted(toolName, lowPassImputationPipelineVersion),
        toolConfig.methodNameWithPipelineVersion());
    assertEquals(
        "%s_v%s".formatted(lowPassPipelineName.getValue(), lowPassImputationPipelineVersion),
        toolConfig.dataTableEntityName());
    assertEquals(pipelineInputDefinitions, toolConfig.inputDefinitions());
    assertEquals(pipelineOutputDefinitions, toolConfig.outputDefinitions());
    assertEquals(useCallCachingPipeline, toolConfig.callCache());
    assertEquals(monitoringScriptPath, toolConfig.monitoringScriptPath());
    assertEquals(deleteIntermediateFilesPipeline, toolConfig.deleteIntermediateOutputFiles());
    assertEquals(memoryRetryMultiplierPipeline, toolConfig.memoryRetryMultiplier());
    assertEquals(pollingIntervalSecondsPipeline, toolConfig.pollingIntervalSeconds());
  }

  @Test
  void testGetPipelineMainToolConfigWithUnsupportedPipeline() {
    // Arrange
    Pipeline pipeline = new Pipeline();

    // Act & Assert
    assertThrows(
        IllegalArgumentException.class,
        () -> toolConfigService.getPipelineMainToolConfig(pipeline),
        "Unsupported pipeline type: null");
  }

  @Test
  void testGetQuotaConsumedToolConfig() {
    Pipeline pipeline = TestUtils.TEST_ARRAY_IMPUTATION_PIPELINE_1; // this has pipeline version 0
    pipeline.setToolVersion(toolVersion);
    pipeline.setPipelineInputDefinitions(pipelineInputDefinitions);

    ToolConfig toolConfig = toolConfigService.getQuotaConsumedToolConfig(pipeline);

    List<PipelineOutputDefinition> expectedOutputDefinitions =
        List.of(
            new PipelineOutputDefinition(
                pipeline.getId(),
                "quotaConsumed",
                "quota_consumed",
                null,
                null,
                PipelineVariableTypesEnum.INTEGER,
                true));
    assertEquals("QuotaConsumed", toolConfig.methodName());
    assertEquals(toolVersion, toolConfig.methodVersion());
    assertEquals("QuotaConsumed_v0", toolConfig.methodNameWithPipelineVersion());
    assertEquals("array_imputation_v0", toolConfig.dataTableEntityName());
    assertEquals(pipelineInputDefinitions, toolConfig.inputDefinitions());
    assertEquals(expectedOutputDefinitions, toolConfig.outputDefinitions());
    assertEquals(useCallCachingQuota, toolConfig.callCache());
    assertEquals(monitoringScriptPath, toolConfig.monitoringScriptPath());
    assertTrue(toolConfig.deleteIntermediateOutputFiles());
    assertNull(toolConfig.memoryRetryMultiplier());
    assertEquals(pollingIntervalSecondsQuota, toolConfig.pollingIntervalSeconds());
  }

  @Test
  void getInputQcToolConfig() {
    Pipeline pipeline = TestUtils.TEST_ARRAY_IMPUTATION_PIPELINE_1; // this has pipeline version 0
    pipeline.setToolVersion(toolVersion);
    pipeline.setPipelineInputDefinitions(pipelineInputDefinitions);

    ToolConfig toolConfig = toolConfigService.getInputQcToolConfig(pipeline);

    List<PipelineOutputDefinition> expectedOutputDefinitions =
        List.of(
            new PipelineOutputDefinition(
                pipeline.getId(),
                "passesQc",
                "passes_qc",
                null,
                null,
                PipelineVariableTypesEnum.BOOLEAN,
                true),
            new PipelineOutputDefinition(
                pipeline.getId(),
                "qcMessages",
                "qc_messages",
                null,
                null,
                PipelineVariableTypesEnum.STRING,
                false));
    assertEquals("InputQC", toolConfig.methodName());
    assertEquals(toolVersion, toolConfig.methodVersion());
    assertEquals("InputQC_v0", toolConfig.methodNameWithPipelineVersion());
    assertEquals("array_imputation_v0", toolConfig.dataTableEntityName());
    assertEquals(pipelineInputDefinitions, toolConfig.inputDefinitions());
    assertEquals(expectedOutputDefinitions, toolConfig.outputDefinitions());
    assertEquals(useCallCachingInputQc, toolConfig.callCache());
    assertEquals(monitoringScriptPath, toolConfig.monitoringScriptPath());
    assertTrue(toolConfig.deleteIntermediateOutputFiles());
    assertNull(toolConfig.memoryRetryMultiplier());
    assertEquals(pollingIntervalSecondsInputQc, toolConfig.pollingIntervalSeconds());
  }
}
