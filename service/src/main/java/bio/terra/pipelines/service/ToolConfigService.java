package bio.terra.pipelines.service;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.app.configuration.internal.PipelinesCommonConfiguration;
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

  private final ImputationConfiguration imputationConfiguration;
  private final PipelinesCommonConfiguration pipelinesCommonConfiguration;

  @Autowired
  public ToolConfigService(
      ImputationConfiguration imputationConfiguration,
      PipelinesCommonConfiguration pipelinesCommonConfiguration) {
    this.imputationConfiguration = imputationConfiguration;
    this.pipelinesCommonConfiguration = pipelinesCommonConfiguration;
  }

  public ToolConfig getPipelineMainToolConfig(Pipeline pipeline) {
    // for now we're hard coding the imputationConfiguration here since it's the only pipeline
    if (PipelinesEnum.ARRAY_IMPUTATION.equals(pipeline.getName())) {
      return new ToolConfig(
          pipeline.getToolName(),
          pipeline.getToolVersion(),
          pipeline.getPipelineInputDefinitions(),
          pipeline.getPipelineOutputDefinitions(),
          imputationConfiguration.isUseCallCaching(),
          imputationConfiguration.isDeleteIntermediateFiles(),
          imputationConfiguration.isUseReferenceDisk(),
          imputationConfiguration.getMemoryRetryMultiplier(),
          imputationConfiguration.getCromwellSubmissionPollingIntervalInSeconds());
    }
    throw new IllegalArgumentException("Unsupported pipeline type: " + pipeline.getName());
  }

  public ToolConfig getPipelineInputQcToolConfig(Pipeline pipeline) {
    return new ToolConfig(
        "InputQC",
        pipeline.getToolVersion(),
        pipeline.getPipelineInputDefinitions(),
        List.of(
            new PipelineOutputDefinition(
                null, "passesQc", "passes_qc", PipelineVariableTypesEnum.STRING),
            new PipelineOutputDefinition(
                null, "QcMessages", "qc_messages", PipelineVariableTypesEnum.STRING)),
        pipelinesCommonConfiguration.isInputQcUseCallCaching(),
        true,
        true,
        pipelinesCommonConfiguration.getMemoryRetryMultiplier(),
        pipelinesCommonConfiguration.getInputQcPollingIntervalSeconds());
  }

  public ToolConfig getQuotaConsumedToolConfig(Pipeline pipeline) {
    return new ToolConfig(
        "QuotaConsumed",
        pipeline.getToolVersion(),
        pipeline.getPipelineInputDefinitions(),
        List.of(
            new PipelineOutputDefinition(
                null, "quotaConsumed", "quota_consumed", PipelineVariableTypesEnum.INTEGER)),
        pipelinesCommonConfiguration.isQuotaConsumedUseCallCaching(),
        true,
        false,
        pipelinesCommonConfiguration.getMemoryRetryMultiplier(),
        pipelinesCommonConfiguration.getQuotaConsumedPollingIntervalSeconds());
  }
}
