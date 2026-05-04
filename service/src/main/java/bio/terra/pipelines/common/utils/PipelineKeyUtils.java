package bio.terra.pipelines.common.utils;

/** Utility methods for canonical pipeline keys in the format {pipelineName}_v{version}. */
public final class PipelineKeyUtils {
  private PipelineKeyUtils() {
    throw new IllegalStateException("Attempted to instantiate utility class PipelineKeyUtils");
  }

  public static String buildPipelineKey(PipelinesEnum pipelineName, Integer pipelineVersion) {
    return "%s_v%s".formatted(pipelineName.getValue(), pipelineVersion);
  }

  public static ParsedPipelineKey parsePipelineKey(String pipelineKey) {
    if (pipelineKey == null || pipelineKey.isBlank()) {
      throw new IllegalArgumentException("pipelineKey must be provided");
    }

    int separatorIndex = pipelineKey.lastIndexOf("_v");
    if (separatorIndex <= 0 || separatorIndex == pipelineKey.length() - 2) {
      throw new IllegalArgumentException(
          "pipelineKey must follow format {pipelineName}_v{pipelineVersion}");
    }

    String pipelineNameString = pipelineKey.substring(0, separatorIndex);
    String pipelineVersionString = pipelineKey.substring(separatorIndex + 2);

    PipelinesEnum pipelineName =
        java.util.Arrays.stream(PipelinesEnum.values())
            .filter(value -> value.getValue().equals(pipelineNameString))
            .findFirst()
            .orElseThrow(
                () ->
                    new IllegalArgumentException(
                        "Unknown pipeline name in pipelineKey %s".formatted(pipelineKey)));

    int pipelineVersion;
    try {
      pipelineVersion = Integer.parseInt(pipelineVersionString);
    } catch (NumberFormatException e) {
      throw new IllegalArgumentException(
          "Invalid pipeline version in pipelineKey %s".formatted(pipelineKey), e);
    }

    return new ParsedPipelineKey(pipelineName, pipelineVersion);
  }

  public record ParsedPipelineKey(PipelinesEnum pipelineName, Integer version) {}
}
