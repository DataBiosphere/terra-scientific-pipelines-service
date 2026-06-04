package bio.terra.pipelines.service.pipeline;

import static bio.terra.pipelines.common.utils.PipelineKeyUtils.buildPipelineKey;

import bio.terra.common.exception.NotFoundException;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations.PipelineInputDefinitionConfiguration;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations.PipelineOutputDefinitionConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.model.PipelineDefinition;
import bio.terra.pipelines.model.PipelineInputDefinition;
import bio.terra.pipelines.model.PipelineOutputDefinition;
import jakarta.annotation.PostConstruct;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

/**
 * Service component that provides access to pipeline definitions from the YAML configuration. This
 * component transforms the Spring-bound configuration objects (WdlBasedPipelineConfiguration) into
 * clean domain models (PipelineDefinition).
 *
 * <p>The provider does not manage persistence; it is purely a bridge between the Spring
 * configuration layer and domain models used by business logic.
 */
@Component
public class PipelineDefinitionProvider {
  private static final Logger logger = LoggerFactory.getLogger(PipelineDefinitionProvider.class);

  private final PipelineConfigurations pipelineConfigurations;
  // Keyed by canonical pipeline key (e.g. "array_imputation_v1"). Populated once at startup since
  // YAML config is immutable after Spring binds it.
  private final Map<String, PipelineDefinition> definitionCache = new HashMap<>();

  @Autowired
  public PipelineDefinitionProvider(PipelineConfigurations pipelineConfigurations) {
    this.pipelineConfigurations = pipelineConfigurations;
  }

  @PostConstruct
  void initCache() {
    for (PipelinesEnum name : PipelinesEnum.values()) {
      Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration> versionedConfigs =
          pipelineConfigurations.getPipelines().get(name.getConfigKeyValue());
      if (versionedConfigs == null) continue;
      for (Map.Entry<String, PipelineConfigurations.WdlBasedPipelineConfiguration> entry :
          versionedConfigs.entrySet()) {
        int version = Integer.parseInt(entry.getKey());
        String pipelineKey = buildPipelineKey(name, version);
        definitionCache.put(
            pipelineKey, transformToDomainModel(entry.getValue(), name, version, pipelineKey));
      }
    }
    logger.debug(
        "PipelineDefinitionProvider cache initialized with {} entries", definitionCache.size());
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
    PipelineDefinition definition = definitionCache.get(pipelineKey);
    if (definition == null) {
      throw new NotFoundException(
          "Pipeline definition not found for key '%s'".formatted(pipelineKey));
    }
    return definition;
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

    Map<String, PipelineConfigurations.WdlBasedPipelineConfiguration> versionedConfigs =
        pipelineConfigurations.getPipelines().get(name.getConfigKeyValue());

    if (versionedConfigs == null || versionedConfigs.isEmpty()) {
      return Collections.emptyList();
    }

    return versionedConfigs.keySet().stream()
        .map(Integer::parseInt)
        .sorted(Collections.reverseOrder())
        .map(version -> definitionCache.get(buildPipelineKey(name, version)))
        .filter(Objects::nonNull)
        .toList();
  }

  /**
   * Transform a WdlBasedPipelineConfiguration into a PipelineDefinition domain model.
   *
   * @param configuration the configuration object from the YAML
   * @param name the pipeline name enum
   * @param version the pipeline version
   * @param pipelineKey the canonical pipeline key
   * @return the transformed domain model
   */
  private PipelineDefinition transformToDomainModel(
      PipelineConfigurations.WdlBasedPipelineConfiguration configuration,
      PipelinesEnum name,
      Integer version,
      String pipelineKey) {

    var inputDefinitions = transformInputDefinitions(configuration.getInputDefinitionConfigs());
    var outputDefinitions = transformOutputDefinitions(configuration.getOutputDefinitionConfigs());

    return PipelineDefinition.builder()
        .name(name)
        .version(version)
        .pipelineKey(pipelineKey)
        .displayName(configuration.getDisplayName())
        .description(configuration.getDescription())
        .pipelineType(configuration.getPipelineType())
        .toolName(configuration.getToolName())
        .inputDefinitions(inputDefinitions)
        .outputDefinitions(outputDefinitions)
        .memoryRetryMultiplier(configuration.getMemoryRetryMultiplier())
        .build();
  }

  /**
   * Transform input definition configs into domain models.
   *
   * @param inputConfigs the configuration objects
   * @return list of domain model input definitions
   */
  private List<PipelineInputDefinition> transformInputDefinitions(
      List<PipelineInputDefinitionConfiguration> inputConfigs) {
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
  private List<PipelineOutputDefinition> transformOutputDefinitions(
      List<PipelineOutputDefinitionConfiguration> outputConfigs) {
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
}
