package bio.terra.pipelines.service;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.app.configuration.internal.PipelineCatalogConfigurations;
import bio.terra.pipelines.common.utils.PipelineKeyUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineQuota;
import bio.terra.pipelines.db.entities.PipelineRuntimeMetadata;
import bio.terra.pipelines.model.Pipeline;
import bio.terra.pipelines.model.PipelineInputDefinition;
import bio.terra.pipelines.model.PipelineOutputDefinition;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Stream;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Service;

/** Service for reading pipeline metadata from YAML-backed catalog configuration. */
@Service
public class PipelineCatalogService implements PipelineCatalogRepository {
  private static final Logger logger = LoggerFactory.getLogger(PipelineCatalogService.class);

  private final PipelineCatalogConfigurations pipelineCatalogConfigurations;

  public PipelineCatalogService(PipelineCatalogConfigurations pipelineCatalogConfigurations) {
    this.pipelineCatalogConfigurations = pipelineCatalogConfigurations;
  }

  @Override
  public Optional<PipelineCatalogConfigurations.PipelineDefinition> getDefinition(
      PipelinesEnum pipelineName, Integer version) {
    Map<String, Map<String, PipelineCatalogConfigurations.PipelineDefinition>> definitions =
        pipelineCatalogConfigurations.getDefinitions();
    if (definitions == null) {
      return Optional.empty();
    }

    Map<String, PipelineCatalogConfigurations.PipelineDefinition> versionDefinitions =
        definitions.get(pipelineName.getValue());
    if (versionDefinitions == null) {
      return Optional.empty();
    }

    return Optional.ofNullable(versionDefinitions.get(version.toString()));
  }

  @Override
  public List<ConfiguredPipelineVersion> getConfiguredPipelineVersions() {
    Map<String, Map<String, PipelineCatalogConfigurations.PipelineDefinition>> definitions =
        pipelineCatalogConfigurations.getDefinitions();

    return definitions.entrySet().stream()
        .flatMap(
            pipelineEntry ->
                pipelineEntry.getValue().entrySet().stream()
                    .map(
                        versionEntry ->
                            new ConfiguredPipelineVersion(
                                parsePipelineName(pipelineEntry.getKey()),
                                parseVersion(
                                    versionEntry.getKey(),
                                    parsePipelineName(pipelineEntry.getKey())),
                                versionEntry.getValue())))
        .sorted(
            Comparator.comparing(
                    (ConfiguredPipelineVersion configuredPipelineVersion) ->
                        configuredPipelineVersion.pipelineName().getValue())
                .thenComparing(ConfiguredPipelineVersion::version, Comparator.reverseOrder()))
        .toList();
  }

  @Override
  public Pipeline hydratePipeline(PipelineRuntimeMetadata dbPipeline) {
    Optional<PipelineCatalogConfigurations.PipelineDefinition> definitionOptional =
        getDefinition(dbPipeline.getName(), dbPipeline.getVersion());

    if (definitionOptional.isEmpty()) {
      logger.warn(
          "No YAML catalog definition found for pipeline {} version {}; falling back to DB metadata",
          dbPipeline.getName(),
          dbPipeline.getVersion());
      return new Pipeline(
          dbPipeline.getName(),
          dbPipeline.getVersion(),
          dbPipeline.isHidden(),
          null,
          null,
          null,
          null,
          dbPipeline.getToolVersion(),
          dbPipeline.getWorkspaceBillingProject(),
          dbPipeline.getWorkspaceName(),
          dbPipeline.getWorkspaceStorageContainerName(),
          dbPipeline.getWorkspaceGoogleProject(),
          List.of(),
          List.of());
    }

    PipelineCatalogConfigurations.PipelineDefinition definition = definitionOptional.get();

    Pipeline hydratedPipeline =
        new Pipeline(
            dbPipeline.getName(),
            dbPipeline.getVersion(),
            dbPipeline.isHidden(),
            definition.getDisplayName(),
            definition.getDescription(),
            definition.getPipelineType(),
            definition.getToolName(),
            dbPipeline.getToolVersion(),
            dbPipeline.getWorkspaceBillingProject(),
            dbPipeline.getWorkspaceName(),
            dbPipeline.getWorkspaceStorageContainerName(),
            dbPipeline.getWorkspaceGoogleProject(),
            buildInputDefinitions(
                PipelineKeyUtils.buildPipelineKey(dbPipeline.getName(), dbPipeline.getVersion()),
                definition.getInputDefinitions()),
            buildOutputDefinitions(
                PipelineKeyUtils.buildPipelineKey(dbPipeline.getName(), dbPipeline.getVersion()),
                definition.getOutputDefinitions()));
    return hydratedPipeline;
  }

  @Override
  public Pipeline buildPipeline(
      PipelinesEnum pipelineName, Integer version, PipelineRuntimeMetadata dbPipelineOrNull) {
    PipelineCatalogConfigurations.PipelineDefinition definition =
        getDefinition(pipelineName, version)
            .orElseThrow(
                () ->
                    new InternalServerErrorException(
                        "No YAML catalog definition found for pipeline %s version %s"
                            .formatted(pipelineName, version)));

    Pipeline pipeline =
        new Pipeline(
            pipelineName,
            version,
            dbPipelineOrNull != null && dbPipelineOrNull.isHidden(),
            definition.getDisplayName(),
            definition.getDescription(),
            definition.getPipelineType(),
            definition.getToolName(),
            dbPipelineOrNull == null ? null : dbPipelineOrNull.getToolVersion(),
            dbPipelineOrNull == null ? null : dbPipelineOrNull.getWorkspaceBillingProject(),
            dbPipelineOrNull == null ? null : dbPipelineOrNull.getWorkspaceName(),
            dbPipelineOrNull == null ? null : dbPipelineOrNull.getWorkspaceStorageContainerName(),
            dbPipelineOrNull == null ? null : dbPipelineOrNull.getWorkspaceGoogleProject(),
            buildInputDefinitions(
                PipelineKeyUtils.buildPipelineKey(pipelineName, version),
                definition.getInputDefinitions()),
            buildOutputDefinitions(
                PipelineKeyUtils.buildPipelineKey(pipelineName, version),
                definition.getOutputDefinitions()));

    return pipeline;
  }

  @Override
  public PipelineQuota getPipelineQuota(PipelinesEnum pipelineName) {
    PipelineCatalogConfigurations.QuotaDefinition quotaDefinition =
        getLatestDefinition(pipelineName).getQuota();
    return new PipelineQuota(
        pipelineName,
        quotaDefinition.getDefaultQuota(),
        quotaDefinition.getMinQuotaConsumed(),
        quotaDefinition.getQuotaUnits());
  }

  private PipelineCatalogConfigurations.PipelineDefinition getLatestDefinition(
      PipelinesEnum pipelineName) {
    Map<String, PipelineCatalogConfigurations.PipelineDefinition> definitionsByVersion =
        pipelineCatalogConfigurations.getDefinitions().get(pipelineName.getValue());
    if (definitionsByVersion == null || definitionsByVersion.isEmpty()) {
      throw new InternalServerErrorException(
          "No YAML catalog definitions found for pipeline %s".formatted(pipelineName));
    }

    return definitionsByVersion.entrySet().stream()
        .max(Comparator.comparingInt(entry -> parseVersion(entry.getKey(), pipelineName)))
        .orElseThrow(
            () ->
                new InternalServerErrorException(
                    "No YAML catalog definitions found for pipeline %s".formatted(pipelineName)))
        .getValue();
  }

  private int parseVersion(String version, PipelinesEnum pipelineName) {
    try {
      return Integer.parseInt(version);
    } catch (NumberFormatException e) {
      throw new InternalServerErrorException(
          "Invalid pipeline version %s configured for pipeline %s".formatted(version, pipelineName),
          e);
    }
  }

  private PipelinesEnum parsePipelineName(String pipelineName) {
    return Stream.of(PipelinesEnum.values())
        .filter(value -> value.getValue().equals(pipelineName))
        .findFirst()
        .orElseThrow(
            () ->
                new InternalServerErrorException(
                    "Invalid pipeline name %s configured in YAML catalog".formatted(pipelineName)));
  }

  private List<PipelineInputDefinition> buildInputDefinitions(
      String pipelineKey, List<PipelineCatalogConfigurations.InputDefinition> inputDefinitions) {
    return inputDefinitions.stream()
        .map(
            input ->
                new PipelineInputDefinition(
                    pipelineKey,
                    input.getName(),
                    input.getWdlVariableName(),
                    input.getDisplayName(),
                    input.getDescription(),
                    input.getType(),
                    input.getFileSuffix(),
                    input.getRequired(),
                    input.getUserProvided(),
                    input.getDefaultValue(),
                    input.getMinValue(),
                    input.getMaxValue()))
        .toList();
  }

  private List<PipelineOutputDefinition> buildOutputDefinitions(
      String pipelineKey, List<PipelineCatalogConfigurations.OutputDefinition> outputDefinitions) {
    return outputDefinitions.stream()
        .map(
            output ->
                new PipelineOutputDefinition(
                    pipelineKey,
                    output.getName(),
                    output.getWdlVariableName(),
                    output.getDisplayName(),
                    output.getDescription(),
                    output.getType(),
                    output.getRequired()))
        .toList();
  }

  public record ConfiguredPipelineVersion(
      PipelinesEnum pipelineName,
      Integer version,
      PipelineCatalogConfigurations.PipelineDefinition definition) {}
}
