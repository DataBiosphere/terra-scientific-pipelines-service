package bio.terra.pipelines.stairway.steps.exception;

import bio.terra.common.exception.InternalServerErrorException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class InternalStepException extends InternalServerErrorException {

  public InternalStepException(String errorMsg) {
    super(
        "An error occurred while running the job. %s Please contact support for help."
            .formatted(errorMsg));
  }
}
