package bio.terra.pipelines.app.configuration.internal;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import jakarta.annotation.PostConstruct;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import lombok.Getter;
import lombok.Setter;
import org.springframework.boot.context.properties.ConfigurationProperties;

/**
 * Configuration-backed pipeline catalog metadata.
 *
 * <p>This model is intended to hold pipeline metadata that is being moved from database tables into
 * YAML configuration.
 */
@Getter
@Setter
@ConfigurationProperties(prefix = "pipelines.catalog")
public class PipelineCatalogConfigurations {
  // pipelineName -> version -> definition
  private Map<String, Map<String, PipelineDefinition>> definitions;

  @PostConstruct
  public void validate() {
    if (definitions == null || definitions.isEmpty()) {
      throw new IllegalArgumentException("pipelines.catalog.definitions must not be empty");
    }

    for (Map.Entry<String, Map<String, PipelineDefinition>> pipelineEntry :
        definitions.entrySet()) {
      String pipelineName = pipelineEntry.getKey();
      requireNonBlank(
          pipelineName, "Pipeline name in pipelines.catalog.definitions must not be blank");

      Map<String, PipelineDefinition> versions = pipelineEntry.getValue();
      if (versions == null || versions.isEmpty()) {
        throw new IllegalArgumentException(
            "Pipeline %s must define at least one version".formatted(pipelineName));
      }

      for (Map.Entry<String, PipelineDefinition> versionEntry : versions.entrySet()) {
        String version = versionEntry.getKey();
        requireNonBlank(
            version,
            "Pipeline version in pipelines.catalog.definitions[%s] must not be blank"
                .formatted(pipelineName));

        PipelineDefinition definition = versionEntry.getValue();
        if (definition == null) {
          throw new IllegalArgumentException(
              "Pipeline definition must be present for %s version %s"
                  .formatted(pipelineName, version));
        }

        definition.validate(pipelineName, version);
      }
    }
  }

  @Getter
  @Setter
  public static class PipelineDefinition {
    private String displayName;
    private String description;
    private String pipelineType;
    private String toolName;
    private QuotaDefinition quota;
    private List<InputDefinition> inputDefinitions;
    private List<OutputDefinition> outputDefinitions;

    private void validate(String pipelineName, String version) {
      requireNonBlank(
          displayName,
          "displayName is required for %s version %s".formatted(pipelineName, version));
      requireNonBlank(
          description,
          "description is required for %s version %s".formatted(pipelineName, version));
      requireNonBlank(
          pipelineType,
          "pipelineType is required for %s version %s".formatted(pipelineName, version));
      requireNonBlank(
          toolName, "toolName is required for %s version %s".formatted(pipelineName, version));

      if (quota == null) {
        throw new IllegalArgumentException(
            "quota is required for %s version %s".formatted(pipelineName, version));
      }
      quota.validate(pipelineName, version);

      if (inputDefinitions == null || inputDefinitions.isEmpty()) {
        throw new IllegalArgumentException(
            "inputDefinitions must not be empty for %s version %s"
                .formatted(pipelineName, version));
      }
      if (outputDefinitions == null || outputDefinitions.isEmpty()) {
        throw new IllegalArgumentException(
            "outputDefinitions must not be empty for %s version %s"
                .formatted(pipelineName, version));
      }

      validateUniqueNames(inputDefinitions, "input", pipelineName, version);
      validateUniqueNames(outputDefinitions, "output", pipelineName, version);

      inputDefinitions.forEach(input -> input.validate(pipelineName, version));
      outputDefinitions.forEach(output -> output.validate(pipelineName, version));
    }
  }

  @Getter
  @Setter
  public static class QuotaDefinition {
    private Integer defaultQuota;
    private Integer minQuotaConsumed;
    private QuotaUnitsEnum quotaUnits;

    private void validate(String pipelineName, String version) {
      requireNonNull(
          defaultQuota,
          "quota.defaultQuota is required for %s version %s".formatted(pipelineName, version));
      requireNonNull(
          minQuotaConsumed,
          "quota.minQuotaConsumed is required for %s version %s".formatted(pipelineName, version));
      requireNonNull(
          quotaUnits,
          "quota.quotaUnits is required for %s version %s".formatted(pipelineName, version));

      if (defaultQuota < 0) {
        throw new IllegalArgumentException(
            "quota.defaultQuota must be >= 0 for %s version %s".formatted(pipelineName, version));
      }
      if (minQuotaConsumed < 0) {
        throw new IllegalArgumentException(
            "quota.minQuotaConsumed must be >= 0 for %s version %s"
                .formatted(pipelineName, version));
      }
    }
  }

  @Getter
  @Setter
  public abstract static class VariableDefinition {
    private String name;
    private String wdlVariableName;
    private String displayName;
    private String description;
    private PipelineVariableTypesEnum type;
    private Boolean required;

    protected void validateShared(String variableKind, String pipelineName, String version) {
      requireNonBlank(
          name,
          "%s name is required for %s version %s".formatted(variableKind, pipelineName, version));
      requireNonBlank(
          wdlVariableName,
          "%s wdlVariableName is required for %s version %s"
              .formatted(variableKind, pipelineName, version));
      requireNonNull(
          type,
          "%s type is required for %s version %s".formatted(variableKind, pipelineName, version));
      requireNonNull(
          required,
          "%s required is required for %s version %s"
              .formatted(variableKind, pipelineName, version));
    }
  }

  @Getter
  @Setter
  public static class InputDefinition extends VariableDefinition {
    private String fileSuffix;
    private Boolean userProvided;
    private Boolean expectsCustomValue;
    private String defaultValue;
    private Double minValue;
    private Double maxValue;

    private void validate(String pipelineName, String version) {
      validateShared("input", pipelineName, version);
      requireNonNull(
          userProvided,
          "input userProvided is required for %s version %s".formatted(pipelineName, version));
      requireNonNull(
          expectsCustomValue,
          "input expectsCustomValue is required for %s version %s"
              .formatted(pipelineName, version));

      if (userProvided && expectsCustomValue) {
        throw new IllegalArgumentException(
            "input %s cannot be both userProvided and expectsCustomValue in %s version %s"
                .formatted(getName(), pipelineName, version));
      }

      if ((getType() == PipelineVariableTypesEnum.FILE
              || getType() == PipelineVariableTypesEnum.FILE_ARRAY
              || getType() == PipelineVariableTypesEnum.MANIFEST)
          && (fileSuffix == null || fileSuffix.isBlank())) {
        throw new IllegalArgumentException(
            "input %s must define fileSuffix for file-like types in %s version %s"
                .formatted(getName(), pipelineName, version));
      }

      if (minValue != null && maxValue != null && minValue > maxValue) {
        throw new IllegalArgumentException(
            "input %s has minValue > maxValue in %s version %s"
                .formatted(getName(), pipelineName, version));
      }
    }
  }

  @Getter
  @Setter
  public static class OutputDefinition extends VariableDefinition {
    private void validate(String pipelineName, String version) {
      validateShared("output", pipelineName, version);
    }
  }

  private static void validateUniqueNames(
      List<? extends VariableDefinition> definitions,
      String variableKind,
      String pipelineName,
      String version) {
    Set<String> names = new HashSet<>();
    Set<String> wdlVariableNames = new HashSet<>();

    for (VariableDefinition definition : definitions) {
      if (!names.add(definition.getName())) {
        throw new IllegalArgumentException(
            "Duplicate %s name %s for %s version %s"
                .formatted(variableKind, definition.getName(), pipelineName, version));
      }
      if (!wdlVariableNames.add(definition.getWdlVariableName())) {
        throw new IllegalArgumentException(
            "Duplicate %s wdlVariableName %s for %s version %s"
                .formatted(variableKind, definition.getWdlVariableName(), pipelineName, version));
      }
    }
  }

  private static void requireNonBlank(String value, String message) {
    if (value == null || value.isBlank()) {
      throw new IllegalArgumentException(message);
    }
  }

  private static void requireNonNull(Object value, String message) {
    if (value == null) {
      throw new IllegalArgumentException(message);
    }
  }
}
