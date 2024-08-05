package bio.terra.pipelines.dependencies.rawls;

import bio.terra.rawls.client.ApiException;

public class RawlsServiceApiException extends RawlsServiceException {

  public RawlsServiceApiException(ApiException exception) {
    super("Rawls returned an unsuccessful status code", exception);
  }

  public RawlsServiceApiException(String message) {
    super(message);
  }
}
