package bio.terra.pipelines.stairway.steps.exception;

import bio.terra.common.exception.InternalServerErrorException;
import java.util.List;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class InternalStepException extends InternalServerErrorException {

  public InternalStepException(String errorMsg) {
    super(
        "An error occurred while running the job. %s Please contact support for help."
            .formatted(errorMsg));
  }

  /**
   * Deserialization constructor used by {@link
   * bio.terra.pipelines.dependencies.stairway.StairwayExceptionSerializer}. The message is already
   * fully formatted, so it is passed through to the superclass without re-wrapping.
   */
  public InternalStepException(String message, List<String> causes) {
    super(message, causes);
  }
}
