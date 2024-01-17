package bio.terra.pipelines.dependencies.wds;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.databiosphere.workspacedata.api.GeneralWdsInformationApi;
import org.databiosphere.workspacedata.api.RecordsApi;
import org.databiosphere.workspacedata.api.SchemaApi;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class WdsClientTest extends BaseEmbeddedDbTest {
  @Autowired WdsClient wdsClient;

  String wdsBaseUri = "wdsBaseUri";
  String auth_token = "authToken";

  @Test
  void TestWdsClientApis() {

    RecordsApi recordsApi = wdsClient.recordsApi(wdsBaseUri, auth_token);

    assertEquals(wdsBaseUri, recordsApi.getApiClient().getBasePath());
    assertTrue(recordsApi.getApiClient().isDebugging());

    GeneralWdsInformationApi generalWdsInformationApi =
        wdsClient.generalWdsInformationApi(wdsBaseUri, auth_token);

    assertEquals(wdsBaseUri, generalWdsInformationApi.getApiClient().getBasePath());
    assertTrue(generalWdsInformationApi.getApiClient().isDebugging());

    SchemaApi schemaApi = wdsClient.schemaApi(wdsBaseUri, auth_token);

    assertEquals(wdsBaseUri, schemaApi.getApiClient().getBasePath());
    assertTrue(schemaApi.getApiClient().isDebugging());
  }
}
