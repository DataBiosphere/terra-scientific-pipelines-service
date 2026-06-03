package bio.terra.pipelines.common.utils;

import bio.terra.common.exception.NotFoundException;
import java.util.Arrays;

public class PipelineKeyUtils {

  // ---------------------------------------------------------------------------
  // Pipeline-key utilities  (canonical format: {pipeline_name}_v{version})
  // ---------------------------------------------------------------------------

  private static final String PIPELINE_NOT_FOUND_MESSAGE = "Pipeline not found for pipelineKey %s";

  /**
   * Build the canonical pipeline key from a name enum and integer version.
   *
   * @param name the pipeline name enum
   * @param version the pipeline version
   * @return key in the form {@code array_imputation_v2}
   */
  public static String buildPipelineKey(PipelinesEnum name, int version) {
    return "%s_v%d".formatted(name.getLowerCaseValue(), version);
  }

  /**
   * Extract the {@link PipelinesEnum} from a canonical pipeline key.
   *
   * @param pipelineKey e.g. {@code array_imputation_v2}
   * @return the matching enum constant
   * @throws NotFoundException if the key cannot be parsed or the name is unrecognised
   */
  public static PipelinesEnum enumFromPipelineKey(String pipelineKey) {
    String nameValue = nameValueFromPipelineKey(pipelineKey);
    return Arrays.stream(PipelinesEnum.values())
        .filter(p -> p.getLowerCaseValue().equals(nameValue))
        .findFirst()
        .orElseThrow(
            () -> new NotFoundException(PIPELINE_NOT_FOUND_MESSAGE.formatted(pipelineKey)));
  }

  /**
   * Extract the integer version from a canonical pipeline key.
   *
   * @param pipelineKey e.g. {@code array_imputation_v2}
   * @return the version number
   * @throws NotFoundException if the key cannot be parsed
   */
  public static int versionFromPipelineKey(String pipelineKey) {
    int separatorIndex = lastSeparatorIndex(pipelineKey);
    try {
      return Integer.parseInt(pipelineKey.substring(separatorIndex + 2));
    } catch (NumberFormatException e) {
      throw new NotFoundException(PIPELINE_NOT_FOUND_MESSAGE.formatted(pipelineKey));
    }
  }

  // -- private helpers --

  private static String nameValueFromPipelineKey(String pipelineKey) {
    int separatorIndex = lastSeparatorIndex(pipelineKey);
    return pipelineKey.substring(0, separatorIndex);
  }

  private static int lastSeparatorIndex(String pipelineKey) {
    if (pipelineKey == null) {
      throw new NotFoundException(PIPELINE_NOT_FOUND_MESSAGE.formatted("null"));
    }
    int idx = pipelineKey.lastIndexOf("_v");
    if (idx < 1 || idx >= pipelineKey.length() - 2) {
      throw new NotFoundException(PIPELINE_NOT_FOUND_MESSAGE.formatted(pipelineKey));
    }
    return idx;
  }
}
