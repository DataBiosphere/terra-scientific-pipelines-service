package bio.terra.pipelines.dependencies.cbas;

import bio.terra.cbas.client.ApiException;

public class CbasServiceApiException extends CbasServiceException {

  public CbasServiceApiException(ApiException exception) {
    super("Cbas returned an unsuccessful status code", exception);
  }

  public CbasServiceApiException(String message) {
    super(message);
  }
}
