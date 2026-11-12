package bio.terra.pipelines.dependencies.rawls;

import bio.terra.rawls.client.ApiException;

@SuppressWarnings("java:S110") // Disable "Inheritance tree of classes should not be too deep"
public class RawlsServiceApiException extends RawlsServiceException {

  public RawlsServiceApiException(ApiException exception) {
    super("A dependent service returned an unsuccessful status code", exception);
  }

  public RawlsServiceApiException(String message) {
    super(message);
  }
}
