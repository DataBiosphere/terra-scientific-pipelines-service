package bio.terra.pipelines.dependencies.leonardo;

import org.broadinstitute.dsde.workbench.client.leonardo.ApiException;

public class LeonardoServiceApiException extends LeonardoServiceException {

  public LeonardoServiceApiException(ApiException exception) {
    super("Leonardo returned an unsuccessful status code", exception);
  }

  public LeonardoServiceApiException(String message) {
    super("Leonardo returned an unsuccessful status code");
  }
}
