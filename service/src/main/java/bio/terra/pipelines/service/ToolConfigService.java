package bio.terra.pipelines.service;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.model.Pipeline;
import bio.terra.pipelines.model.PipelineOutputDefinition;
import bio.terra.pipelines.stairway.steps.utils.ToolConfig;
import java.util.List;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

/** Service to encapsulate logic used to generate Tool Configs for pipelines */
@Service
public class ToolConfigService {

  private final PipelineConfigurations pipelineConfigurations;

  @Autowired
  public ToolConfigService(PipelineConfigurations pipelineConfigurations) {
    this.pipelineConfigurations = pipelineConfigurations;
  }

  public static final String INPUT_QC_METHOD_NAME = "InputQC";
  public static final String QUOTA_CONSUMED_METHOD_NAME = "QuotaConsumed";

  /** Get the ToolConfig for the main analysis method/workflow for a given pipeline */
  public ToolConfig getPipelineMainToolConfig(Pipeline pipeline) {
    PipelineConfigurations.CommonConfiguration commonConfiguration =
        pipelineConfigurations.getCommon();
    String toolNameWithPipelineVersion =
        appendPipelineVersion(pipeline.getToolName(), pipeline.getVersion());
    PipelineConfigurations.WdlBasedPipelineConfiguration pipelineConfiguration =
        pipelineConfigurations.getPipelineConfiguration(pipeline.getPipelineKey());
    PipelineConfigurations.PipelineMetadataConfiguration metadata =
        pipelineConfiguration.getMetadata();
    return new ToolConfig(
        pipeline.getToolName(),
        pipeline.getToolVersion(),
        toolNameWithPipelineVersion,
        getDataTableEntityNameForToolConfig(pipeline),
        pipeline.getInputDefinitions(),
        pipeline.getOutputDefinitions(),
        commonConfiguration.isMainToolUseCallCaching(),
        commonConfiguration.getMonitoringScriptPath(),
        commonConfiguration.isMainToolDeleteIntermediateFiles(),
        metadata.getMemoryRetryMultiplier(),
        commonConfiguration.getMainToolPollingIntervalSeconds());
  }

  /**
   * Get the ToolConfig for the InputQC method for a given pipeline. InputQC wdls are expected to
   * take the same inputs as the main pipeline wdl, and output two variables: boolean passes_qc and
   * string qc_messages. If passes_qc is true, qc_messages should be an empty string. If passes_qc
   * is false, qc_messages should contain user-facing, actionable messages about the qc failure(s).
   */
  public ToolConfig getInputQcToolConfig(Pipeline pipeline) {
    String methodNameWithPipelineVersion =
        appendPipelineVersion(INPUT_QC_METHOD_NAME, pipeline.getVersion());
    PipelineConfigurations.CommonConfiguration commonConfiguration =
        pipelineConfigurations.getCommon();
    return new ToolConfig(
        INPUT_QC_METHOD_NAME,
        pipeline.getToolVersion(),
        methodNameWithPipelineVersion,
        getDataTableEntityNameForToolConfig(pipeline),
        pipeline.getInputDefinitions(),
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
                .build()),
        commonConfiguration.isInputQcUseCallCaching(),
        commonConfiguration.getMonitoringScriptPath(),
        true,
        null, // no memory retry multiplier
        commonConfiguration.getInputQcPollingIntervalSeconds());
  }

  /**
   * Get the ToolConfig for the QuotaConsumed method for a given pipeline. QuotaConsumed wdls are
   * expected to take the same inputs as the main pipeline wdl, and output one variable: integer
   * quota_consumed.
   */
  public ToolConfig getQuotaConsumedToolConfig(Pipeline pipeline) {
    String methodNameWithPipelineVersion =
        appendPipelineVersion(QUOTA_CONSUMED_METHOD_NAME, pipeline.getVersion());
    PipelineConfigurations.CommonConfiguration commonConfiguration =
        pipelineConfigurations.getCommon();
    return new ToolConfig(
        QUOTA_CONSUMED_METHOD_NAME,
        pipeline.getToolVersion(),
        methodNameWithPipelineVersion,
        getDataTableEntityNameForToolConfig(pipeline),
        pipeline.getInputDefinitions(),
        List.of(
            PipelineOutputDefinition.builder()
                .name("quotaConsumed")
                .wdlVariableName("quota_consumed")
                .type(PipelineVariableTypesEnum.INTEGER)
                .isRequired(true)
                .build()),
        commonConfiguration.isQuotaConsumedUseCallCaching(),
        commonConfiguration.getMonitoringScriptPath(),
        true,
        null, // no memory retry multiplier
        commonConfiguration.getQuotaConsumedPollingIntervalSeconds());
  }

  /**
   * Helper method to construct a {inputString}_v{version} string.
   *
   * @param inputString method name
   * @param pipelineVersion integer
   * @return {inputString}_v{pipelineVersion}
   */
  private String appendPipelineVersion(String inputString, Integer pipelineVersion) {
    return "%s_v%s".formatted(inputString, pipelineVersion);
  }

  /** Helper method to construct the data table entity name for a given pipeline version */
  private String getDataTableEntityNameForToolConfig(Pipeline pipeline) {
    return appendPipelineVersion(pipeline.getName().getLowerCaseValue(), pipeline.getVersion());
  }
}
