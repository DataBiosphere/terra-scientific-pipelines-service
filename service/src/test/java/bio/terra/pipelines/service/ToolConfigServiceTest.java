package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.app.configuration.internal.PipelinesCommonConfiguration;
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
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

class ToolConfigServiceTest extends BaseTest {

  @Mock private ImputationConfiguration imputationConfiguration;
  @Mock private PipelinesCommonConfiguration pipelinesCommonConfiguration;

  private ToolConfigService toolConfigService;

  private final PipelinesEnum pipelineName = PipelinesEnum.ARRAY_IMPUTATION;
  private final String toolName = TestUtils.TEST_TOOL_NAME_1;
  private final String toolVersion = TestUtils.TEST_TOOL_VERSION_1;
  private final List<PipelineInputDefinition> pipelineInputDefinitions =
      TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST;
  private final List<PipelineOutputDefinition> pipelineOutputDefinitions =
      TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST;
  private final boolean useCallCachingPipeline = true;
  private final boolean deleteIntermediateFilesPipeline = false;
  private final boolean useReferenceDiskPipeline = true;
  private final BigDecimal memoryRetryMultiplierPipeline = BigDecimal.valueOf(1.5);
  private final Long pollingIntervalSecondsPipeline = 30L;
  private final Long pollingIntervalSecondsQuota = 5L;
  private final boolean useCallCachingQuota = false;

  @BeforeEach
  void setUp() {
    MockitoAnnotations.openMocks(this);
    toolConfigService =
        new ToolConfigService(imputationConfiguration, pipelinesCommonConfiguration);

    // mock imputationConfig
    when(imputationConfiguration.isUseCallCaching()).thenReturn(useCallCachingPipeline);
    when(imputationConfiguration.isDeleteIntermediateFiles())
        .thenReturn(deleteIntermediateFilesPipeline);
    when(imputationConfiguration.isUseReferenceDisk()).thenReturn(useReferenceDiskPipeline);
    when(imputationConfiguration.getMemoryRetryMultiplier())
        .thenReturn(memoryRetryMultiplierPipeline);
    when(imputationConfiguration.getCromwellSubmissionPollingIntervalInSeconds())
        .thenReturn(pollingIntervalSecondsPipeline);

    // mock pipelinesCommonConfiguration
    when(pipelinesCommonConfiguration.getQuotaConsumedPollingIntervalSeconds())
        .thenReturn(pollingIntervalSecondsQuota);
    when(pipelinesCommonConfiguration.isQuotaConsumedUseCallCaching())
        .thenReturn(useCallCachingQuota);
  }

  @Test
  void testGetPipelineMainToolConfig_WithImputationPipeline() {
    // create pipeline
    Pipeline pipeline = new Pipeline();
    pipeline.setName(pipelineName);
    pipeline.setToolName(toolName);
    pipeline.setToolVersion(toolVersion);
    pipeline.setPipelineInputDefinitions(pipelineInputDefinitions);
    pipeline.setPipelineOutputDefinitions(pipelineOutputDefinitions);

    // create main tool config
    ToolConfig toolConfig = toolConfigService.getPipelineMainToolConfig(pipeline);

    // check values
    assertEquals(toolName, toolConfig.methodName());
    assertEquals(toolVersion, toolConfig.methodVersion());
    assertEquals(pipelineInputDefinitions, toolConfig.inputDefinitions());
    assertEquals(pipelineOutputDefinitions, toolConfig.outputDefinitions());
    assertEquals(useCallCachingPipeline, toolConfig.callCache());
    assertEquals(deleteIntermediateFilesPipeline, toolConfig.deleteIntermediateOutputFiles());
    assertEquals(useReferenceDiskPipeline, toolConfig.useReferenceDisks());
    assertEquals(memoryRetryMultiplierPipeline, toolConfig.memoryRetryMultiplier());
    assertEquals(pollingIntervalSecondsPipeline, toolConfig.pollingIntervalSeconds());
  }

  @Test
  void testGetPipelineMainToolConfig_WithUnsupportedPipeline() {
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
    Pipeline pipeline = new Pipeline();
    pipeline.setToolVersion(toolVersion);
    pipeline.setPipelineInputDefinitions(pipelineInputDefinitions);

    ToolConfig toolConfig = toolConfigService.getQuotaConsumedToolConfig(pipeline);

    List<PipelineOutputDefinition> expectedOutputDefinitions =
        List.of(
            new PipelineOutputDefinition(
                null, "quotaConsumed", "quota_consumed", PipelineVariableTypesEnum.INTEGER));
    assertEquals("QuotaConsumed", toolConfig.methodName());
    assertEquals(toolVersion, toolConfig.methodVersion());
    assertEquals(pipelineInputDefinitions, toolConfig.inputDefinitions());
    assertEquals(expectedOutputDefinitions, toolConfig.outputDefinitions());
    assertEquals(useCallCachingQuota, toolConfig.callCache());
    assertTrue(toolConfig.deleteIntermediateOutputFiles());
    assertFalse(toolConfig.useReferenceDisks());
    assertNull(toolConfig.memoryRetryMultiplier());
    assertEquals(pollingIntervalSecondsQuota, toolConfig.pollingIntervalSeconds());
  }
}
