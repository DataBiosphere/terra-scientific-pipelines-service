package bio.terra.pipelines.service;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
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

  /** Get the ToolConfig for the main analysis method/workflow for a given pipeline */
  public ToolConfig getPipelineMainToolConfig(Pipeline pipeline) {
    // for now we're hard coding the imputationConfiguration here since it's the only pipeline
    if (PipelinesEnum.ARRAY_IMPUTATION.equals(pipeline.getName())) {
      PipelineConfigurations.ArrayImputationConfig arrayImputationConfiguration =
          pipelineConfigurations.getArrayImputation().get(pipeline.getVersion().toString());
      String toolNameWithPipelineVersion =
          appendPipelineVersion(pipeline.getToolName(), pipeline.getVersion());
      return new ToolConfig(
          pipeline.getToolName(),
          pipeline.getToolVersion(),
          toolNameWithPipelineVersion,
          getDataTableEntityNameForToolConfig(pipeline),
          pipeline.getPipelineInputDefinitions(),
          pipeline.getPipelineOutputDefinitions(),
          arrayImputationConfiguration.isUseCallCaching(),
          pipelineConfigurations.getCommon().getMonitoringScriptPath(),
          arrayImputationConfiguration.isDeleteIntermediateFiles(),
          arrayImputationConfiguration.getMemoryRetryMultiplier(),
          arrayImputationConfiguration.getCromwellSubmissionPollingIntervalInSeconds());
    }
    throw new IllegalArgumentException("Unsupported pipeline type: " + pipeline.getName());
  }

  /**
   * Get the ToolConfig for the InputQC method for a given pipeline. InputQC wdls are expected to
   * take the same inputs as the main pipeline wdl, and output two variables: boolean passes_qc and
   * string qc_messages. If passes_qc is true, qc_messages should be an empty string. If passes_qc
   * is false, qc_messages should contain user-facing, actionable messages about the qc failure(s).
   */
  public ToolConfig getInputQcToolConfig(Pipeline pipeline) {
    String inputQcMethodName = "InputQC";
    String methodNameWithPipelineVersion =
        appendPipelineVersion(inputQcMethodName, pipeline.getVersion());
    return new ToolConfig(
        inputQcMethodName,
        pipeline.getToolVersion(),
        methodNameWithPipelineVersion,
        getDataTableEntityNameForToolConfig(pipeline),
        pipeline.getPipelineInputDefinitions(),
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
                false)),
        pipelineConfigurations.getCommon().isInputQcUseCallCaching(),
        pipelineConfigurations.getCommon().getMonitoringScriptPath(),
        true,
        null, // no memory retry multiplier
        pipelineConfigurations.getCommon().getInputQcPollingIntervalSeconds());
  }

  /**
   * Get the ToolConfig for the QuotaConsumed method for a given pipeline. QuotaConsumed wdls are
   * expected to take the same inputs as the main pipeline wdl, and output one variable: integer
   * quota_consumed.
   */
  public ToolConfig getQuotaConsumedToolConfig(Pipeline pipeline) {
    String quotaConsumedMethodName = "QuotaConsumed";
    String methodNameWithPipelineVersion =
        appendPipelineVersion(quotaConsumedMethodName, pipeline.getVersion());
    return new ToolConfig(
        quotaConsumedMethodName,
        pipeline.getToolVersion(),
        methodNameWithPipelineVersion,
        getDataTableEntityNameForToolConfig(pipeline),
        pipeline.getPipelineInputDefinitions(),
        List.of(
            new PipelineOutputDefinition(
                pipeline.getId(),
                "quotaConsumed",
                "quota_consumed",
                null,
                null,
                PipelineVariableTypesEnum.INTEGER,
                true)),
        pipelineConfigurations.getCommon().isQuotaConsumedUseCallCaching(),
        pipelineConfigurations.getCommon().getMonitoringScriptPath(),
        true,
        null, // no memory retry multiplier
        pipelineConfigurations.getCommon().getQuotaConsumedPollingIntervalSeconds());
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
    return appendPipelineVersion(pipeline.getName().getValue(), pipeline.getVersion());
  }
}
