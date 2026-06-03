package bio.terra.pipelines.service.pipeline;

import static bio.terra.pipelines.common.utils.PipelineKeyUtils.buildPipelineKey;

import bio.terra.common.exception.NotFoundException;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations.PipelineInputDefinitionConfig;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations.PipelineOutputDefinitionConfig;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.model.PipelineDefinition;
import bio.terra.pipelines.model.PipelineInputDefinition;
import bio.terra.pipelines.model.PipelineOutputDefinition;
import bio.terra.pipelines.model.PipelineQuota;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

/**
 * Service component that provides access to pipeline definitions from the YAML configuration. This
 * component transforms the Spring-bound configuration objects (PipelineConfiguration) into clean
 * domain models (PipelineDefinition).
 *
 * <p>The provider does not manage persistence; it is purely a bridge between the Spring
 * configuration layer and domain models used by business logic.
 */
@Component
public class PipelineDefinitionProvider {
  private static final Logger logger = LoggerFactory.getLogger(PipelineDefinitionProvider.class);

  private final PipelineConfigurations pipelineConfigurations;
  private final Map<String, PipelineDefinition> definitionCache = new HashMap<>();

  @Autowired
  public PipelineDefinitionProvider(PipelineConfigurations pipelineConfigurations) {
    this.pipelineConfigurations = pipelineConfigurations;
  }

  /**
   * Get a pipeline definition by name and version.
   *
   * @param name the pipeline name
   * @param version the pipeline version
   * @return the pipeline definition
   * @throws NotFoundException if the pipeline is not found in the YAML configuration
   */
  public PipelineDefinition getPipelineDefinition(PipelinesEnum name, Integer version) {
    logger.debug("Getting pipeline definition for {} v{}", name.getLowerCaseValue(), version);
    String pipelineKey = buildPipelineKey(name, version);
    PipelineConfigurations.PipelineConfiguration wdlConfig =
        pipelineConfigurations.getPipelineConfiguration(pipelineKey);
    return transformToDomainModel(wdlConfig, name, version, pipelineKey);
  }

  /**
   * Get all pipeline definitions for a given pipeline name, sorted by version descending (highest
   * first). Does not perform any filtering for hidden, since this only looks at YAML-based
   * definitions.
   *
   * @param name the pipeline name
   * @return list of pipeline definitions for all YAML-defined versions, highest version first
   */
  public List<PipelineDefinition> getDefinitionsForPipeline(PipelinesEnum name) {
    logger.debug("Getting all pipeline definitions for {}", name.getLowerCaseValue());

    Map<String, PipelineConfigurations.PipelineConfiguration> versionedConfigs =
        pipelineConfigurations.getPipelines().get(name.getConfigKeyValue());

    if (versionedConfigs == null || versionedConfigs.isEmpty()) {
      return Collections.emptyList();
    }

    return versionedConfigs.entrySet().stream()
        .sorted(
            (e1, e2) ->
                Integer.compare(Integer.parseInt(e2.getKey()), Integer.parseInt(e1.getKey())))
        .map(
            entry -> {
              int version = Integer.parseInt(entry.getKey());
              String pipelineKey = buildPipelineKey(name, version);
              return transformToDomainModel(entry.getValue(), name, version, pipelineKey);
            })
        .toList();
  }

  /**
   * Transform a PipelineConfiguration into a PipelineDefinition domain model.
   *
   * @param configuration the configuration object from the YAML
   * @param name the pipeline name enum
   * @param version the pipeline version
   * @param pipelineKey the canonical pipeline key
   * @return the transformed domain model
   */
  private PipelineDefinition transformToDomainModel(
      PipelineConfigurations.PipelineConfiguration configuration,
      PipelinesEnum name,
      Integer version,
      String pipelineKey) {

    var metadata = configuration.getMetadata();
    var inputDefinitions = transformInputs(configuration.getInputDefinitionConfigs());
    var outputDefinitions = transformOutputs(configuration.getOutputDefinitionConfigs());

    return PipelineDefinition.builder()
        .name(name)
        .version(version)
        .pipelineKey(pipelineKey)
        .displayName(metadata.getDisplayName())
        .description(metadata.getDescription())
        .pipelineType(metadata.getPipelineType())
        .toolName(metadata.getToolName())
        .inputDefinitions(inputDefinitions)
        .outputDefinitions(outputDefinitions)
        .memoryRetryMultiplier(metadata.getMemoryRetryMultiplier())
        .build();
  }

  /**
   * Transform input definition configs into domain models.
   *
   * @param inputConfigs the configuration objects
   * @return list of domain model input definitions
   */
  private List<PipelineInputDefinition> transformInputs(
      List<PipelineInputDefinitionConfig> inputConfigs) {
    if (inputConfigs == null) {
      return Collections.emptyList();
    }
    return inputConfigs.stream()
        .map(
            config ->
                PipelineInputDefinition.builder()
                    .name(config.getName())
                    .wdlVariableName(config.getWdlVariableName())
                    .displayName(config.getDisplayName())
                    .description(config.getDescription())
                    .type(config.getType())
                    .isRequired(config.getIsRequired() != null && config.getIsRequired())
                    .userProvided(config.getUserProvided() != null && config.getUserProvided())
                    .defaultValue(config.getDefaultValue())
                    .minValue(config.getMinValue())
                    .maxValue(config.getMaxValue())
                    .fileSuffix(config.getFileSuffix())
                    .build())
        .collect(Collectors.toList());
  }

  /**
   * Transform output definition configs into domain models.
   *
   * @param outputConfigs the configuration objects
   * @return list of domain model output definitions
   */
  private List<PipelineOutputDefinition> transformOutputs(
      List<PipelineOutputDefinitionConfig> outputConfigs) {
    if (outputConfigs == null) {
      return Collections.emptyList();
    }
    return outputConfigs.stream()
        .map(
            config ->
                PipelineOutputDefinition.builder()
                    .name(config.getName())
                    .wdlVariableName(config.getWdlVariableName())
                    .displayName(config.getDisplayName())
                    .description(config.getDescription())
                    .type(config.getType())
                    .isRequired(config.getIsRequired() != null && config.getIsRequired())
                    .build())
        .collect(Collectors.toList());
  }

  /**
   * Transform a quota config into a domain model.
   *
   * @param quotaConfig the configuration object
   * @return the domain model quota
   */
  private PipelineQuota transformQuota(PipelineConfigurations.PipelineQuotaConfig quotaConfig) {
    return PipelineQuota.builder()
        .defaultQuota(quotaConfig.getDefaultQuota())
        .minQuotaConsumed(quotaConfig.getMinQuotaConsumed())
        .quotaUnits(quotaConfig.getQuotaUnits())
        .build();
  }
}
