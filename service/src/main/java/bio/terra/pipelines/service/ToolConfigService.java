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
      PipelineConfigurations.ImputationConfig imputationConfiguration =
          pipelineConfigurations.getArrayImputation().get(pipeline.getVersion().toString());
      return new ToolConfig(
          pipeline.getToolName(),
          pipeline.getToolVersion(),
          pipeline.getPipelineInputDefinitions(),
          pipeline.getPipelineOutputDefinitions(),
          imputationConfiguration.isUseCallCaching(),
          pipelineConfigurations.getCommon().getMonitoringScriptPath(),
          imputationConfiguration.isDeleteIntermediateFiles(),
          imputationConfiguration.getMemoryRetryMultiplier(),
          imputationConfiguration.getCromwellSubmissionPollingIntervalInSeconds());
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
    return new ToolConfig(
        "InputQC",
        pipeline.getToolVersion(),
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
    return new ToolConfig(
        "QuotaConsumed",
        pipeline.getToolVersion(),
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
}
