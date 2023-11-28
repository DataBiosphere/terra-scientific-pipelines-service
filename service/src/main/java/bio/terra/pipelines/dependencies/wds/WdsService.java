package bio.terra.pipelines.dependencies.wds;

import bio.terra.pipelines.app.configuration.external.WdsServerConfiguration;
import bio.terra.pipelines.dependencies.common.DependencyNotAvailableException;
import bio.terra.pipelines.dependencies.common.HealthCheckWorkspaceApps;
import java.util.List;
import java.util.Objects;
import org.databiosphere.workspacedata.client.ApiException;
import org.databiosphere.workspacedata.model.RecordRequest;
import org.databiosphere.workspacedata.model.RecordResponse;
import org.databiosphere.workspacedata.model.RecordTypeSchema;
import org.databiosphere.workspacedata.model.StatusResponse;
import org.springframework.retry.support.RetryTemplate;
import org.springframework.stereotype.Service;

/** class to encapsulate interacting with WDS client */
@Service
public class WdsService implements HealthCheckWorkspaceApps {

  private final WdsClient wdsClient;
  private final WdsServerConfiguration wdsServerConfiguration;
  private final RetryTemplate listenerResetRetryTemplate;

  public WdsService(
      WdsClient wdsClient,
      WdsServerConfiguration wdsServerConfiguration,
      RetryTemplate listenerResetRetryTemplate) {
    this.wdsClient = wdsClient;
    this.wdsServerConfiguration = wdsServerConfiguration;
    this.listenerResetRetryTemplate = listenerResetRetryTemplate;
  }

  /**
   * query for a record from WDS table
   *
   * @param wdsBaseUri - URI of WDS app, retrieved from Leonardo
   * @param bearerToken - SA token, retrieved from google default credentials
   * @param recordType - WDS table name
   * @param workspaceId - workspace id WDS app is deployed in
   * @param recordId - id for the record
   * @return - returns Record details
   * @throws WdsServiceException ;
   */
  public RecordResponse getRecord(
      String wdsBaseUri, String bearerToken, String recordType, String workspaceId, String recordId)
      throws WdsServiceException {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () ->
            wdsClient
                .recordsApi(wdsBaseUri, bearerToken)
                .getRecord(workspaceId, wdsServerConfiguration.apiV(), recordType, recordId));
  }

  /**
   * Update a record in a table
   *
   * @param wdsBaseUri - URI of WDS app, retrieved from Leonardo
   * @param bearerToken - SA token, retrieved from google default credentials
   * @param request - request for what fields to update
   * @param workspaceId - workspace id WDS app is deployed in
   * @param type - WDS table name
   * @param id - id of record to update
   * @return - the updated record
   * @throws WdsServiceException ;
   */
  public RecordResponse updateRecord(
      String wdsBaseUri,
      String bearerToken,
      RecordRequest request,
      String workspaceId,
      String type,
      String id)
      throws WdsServiceException {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () ->
            wdsClient
                .recordsApi(wdsBaseUri, bearerToken)
                .updateRecord(request, workspaceId, wdsServerConfiguration.apiV(), type, id));
  }

  /**
   * Query for all the tables and their structure for this WDS app
   *
   * @param wdsBaseUri - URI of WDS app, retrieved from Leonardo
   * @param bearerToken - SA token, retrieved from google default credentials
   * @param workspaceId - workspace id WDS app is deployed in
   * @return - list of all the tables in this WDS app
   * @throws WdsServiceException ;
   */
  public List<RecordTypeSchema> querySchema(
      String wdsBaseUri, String bearerToken, String workspaceId) throws WdsServiceException {
    return executionWithRetryTemplate(
        listenerResetRetryTemplate,
        () ->
            wdsClient
                .schemaApi(wdsBaseUri, bearerToken)
                .describeAllRecordTypes(workspaceId, wdsServerConfiguration.apiV()));
  }

  @Override
  public Result checkHealth(String wdsBaseUri, String accessToken) {
    try {
      StatusResponse result =
          wdsClient.generalWdsInformationApi(wdsBaseUri, accessToken).statusGet();
      return new Result(Objects.equals(result.getStatus(), "UP"), result.toString());
    } catch (ApiException e) {
      return new Result(false, e.getMessage());
    }
  }

  interface WdsAction<T> {
    T execute() throws ApiException, DependencyNotAvailableException;
  }

  @SuppressWarnings("java:S125") // The comment here isn't "commented code"
  static <T> T executionWithRetryTemplate(RetryTemplate retryTemplate, WdsAction<T> action)
      throws WdsServiceException {

    // Why all this song and dance to catch exceptions and map them to almost identical exceptions?
    // Because the RetryTemplate's execute function only allows us to declare one Throwable type.
    // So we have a top-level WdsServiceException that we can catch and handle, and then we have
    // subclasses of that exception representing the types of exception that can be thrown. This
    // way, we can keep well typed exceptions (no "catch (Exception e)") and still make use of the
    // retry framework.
    return retryTemplate.execute(
        context -> {
          try {
            return action.execute();
          } catch (ApiException e) {
            throw new WdsServiceApiException(e);
          } catch (DependencyNotAvailableException e) {
            throw new WdsServiceNotAvailableException(e);
          }
        });
  }
}
