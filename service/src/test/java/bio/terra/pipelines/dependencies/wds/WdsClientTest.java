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
  String authToken = "authToken";

  @Test
  void testWdsClientApis() {

    RecordsApi recordsApi = wdsClient.recordsApi(wdsBaseUri, authToken);

    assertEquals(wdsBaseUri, recordsApi.getApiClient().getBasePath());
    assertTrue(recordsApi.getApiClient().isDebugging());

    GeneralWdsInformationApi generalWdsInformationApi =
        wdsClient.generalWdsInformationApi(wdsBaseUri, authToken);

    assertEquals(wdsBaseUri, generalWdsInformationApi.getApiClient().getBasePath());
    assertTrue(generalWdsInformationApi.getApiClient().isDebugging());

    SchemaApi schemaApi = wdsClient.schemaApi(wdsBaseUri, authToken);

    assertEquals(wdsBaseUri, schemaApi.getApiClient().getBasePath());
    assertTrue(schemaApi.getApiClient().isDebugging());
  }
}
