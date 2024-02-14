package bio.terra.pipelines.dependencies.cbas;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.cbas.api.MethodsApi;
import bio.terra.cbas.api.PublicApi;
import bio.terra.cbas.api.RunSetsApi;
import bio.terra.cbas.api.RunsApi;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class CbasClientTest extends BaseEmbeddedDbTest {
  @Autowired CbasClient cbasClient;
  String cbasBaseUri = "cbasBaseUri";
  String authToken = "authToken";

  @Test
  void testWdsClientApis() {

    PublicApi publicApi = cbasClient.publicApi(cbasBaseUri, authToken);

    assertEquals(cbasBaseUri, publicApi.getApiClient().getBasePath());
    assertTrue(publicApi.getApiClient().isDebugging());

    MethodsApi methodsApi = cbasClient.methodsApi(cbasBaseUri, authToken);

    assertEquals(cbasBaseUri, methodsApi.getApiClient().getBasePath());
    assertTrue(methodsApi.getApiClient().isDebugging());

    RunsApi runsApi = cbasClient.runsApi(cbasBaseUri, authToken);

    assertEquals(cbasBaseUri, runsApi.getApiClient().getBasePath());
    assertTrue(runsApi.getApiClient().isDebugging());

    RunSetsApi runSetsApi = cbasClient.runSetsApi(cbasBaseUri, authToken);

    assertEquals(cbasBaseUri, runSetsApi.getApiClient().getBasePath());
    assertTrue(runSetsApi.getApiClient().isDebugging());
  }
}
