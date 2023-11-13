package bio.terra.pipelines.dependencies.leonardo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import bio.terra.pipelines.app.configuration.external.LeonardoServerConfiguration;
import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import com.google.gson.Gson;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.apache.commons.text.StringSubstitutor;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
import org.junit.jupiter.api.Test;

class AppUtilsTest {
  private final String workspaceId = UUID.randomUUID().toString();

  private final LeonardoServerConfiguration leonardoServerConfiguration =
      new LeonardoServerConfiguration(
          "baseuri",
          List.of("WDS", "CROMWELL"),
          List.of("CROMWELL_RUNNER_APP", "CROMWELL"),
          0,
          false);

  private final ImputationConfiguration imputationConfiguration =
      new ImputationConfiguration(workspaceId);
  private final ListAppResponse separatedWdsAppWdsAppName;
  // does not have "wds-" app name and does not have anticipated wds url
  private final ListAppResponse separatedWdsAppNotWdsAppName;
  private final ListAppResponse cromwellListAppResponse;
  private final ListAppResponse combinedWdsInCromwellApp;
  private final ListAppResponse otherNamedCromwellAppOlder;
  private final ListAppResponse galaxyApp;
  private final ListAppResponse otherNamedCromwellApp;
  private final ListAppResponse otherNamedCromwellAppProvisioning;
  private final ListAppResponse separatedWorkflowsApp;
  private final ListAppResponse cromwellRunnerAppResponse;

  private final AppUtils au = new AppUtils(leonardoServerConfiguration, imputationConfiguration);

  @Test
  void findWdsUrlInCombinedApp() throws Exception {
    List<ListAppResponse> apps = List.of(combinedWdsInCromwellApp);

    AppUtils au = new AppUtils(leonardoServerConfiguration, imputationConfiguration);
    assertEquals(anticipatedWdsUrl("wds"), au.findUrlForWds(apps));
  }

  @Test
  void findCromwellUrlInTerraApp() throws Exception {
    List<ListAppResponse> apps = List.of(cromwellListAppResponse);

    AppUtils au = new AppUtils(leonardoServerConfiguration, imputationConfiguration);
    assertEquals(anticipatedCromwellUrl("terra-app"), au.findUrlForCromwell(apps));
  }

  @Test
  void findWdsUrlInWdsApp() throws Exception {
    List<ListAppResponse> apps = List.of(separatedWdsAppWdsAppName);

    AppUtils au = new AppUtils(leonardoServerConfiguration, imputationConfiguration);
    assertEquals(anticipatedWdsUrl("wds"), au.findUrlForWds(apps));
  }

  @Test
  void findWdsAndCromwellInCombinedAppResponse() throws Exception {
    List<ListAppResponse> apps = List.of(separatedWdsAppWdsAppName, cromwellListAppResponse);

    AppUtils au = new AppUtils(leonardoServerConfiguration, imputationConfiguration);
    assertEquals(anticipatedWdsUrl("wds"), au.findUrlForWds(apps));
    assertEquals(anticipatedCromwellUrl("terra-app"), au.findUrlForCromwell(apps));
  }

  @Test
  void preferSpecificallyNamedAppGoodFirst() throws Exception {
    List<ListAppResponse> apps =
        new java.util.ArrayList<>(List.of(separatedWdsAppWdsAppName, separatedWdsAppNotWdsAppName));
    assertEquals(anticipatedWdsUrl("wds"), au.findUrlForWds(apps));
  }

  @Test
  void preferSpecificallyNamedAppGoodSecond() throws Exception {
    List<ListAppResponse> apps =
        new java.util.ArrayList<>(List.of(separatedWdsAppNotWdsAppName, separatedWdsAppWdsAppName));
    assertEquals(anticipatedWdsUrl("wds"), au.findUrlForWds(apps));
  }

  @Test
  void notChooseGalaxyApp() throws Exception {
    List<ListAppResponse> apps =
        new java.util.ArrayList<>(List.of(otherNamedCromwellApp, galaxyApp));
    // Shuffle to make sure the initial ordering isn't relevant:
    Collections.shuffle(apps);
    assertEquals(anticipatedWdsUrl("app1"), au.findUrlForWds(apps));
  }

  @Test
  void preferNewerCreatedApp() throws Exception {
    List<ListAppResponse> apps =
        new java.util.ArrayList<>(List.of(otherNamedCromwellApp, otherNamedCromwellAppOlder));
    // Shuffle to make sure the initial ordering isn't relevant:
    Collections.shuffle(apps);
    assertEquals(anticipatedWdsUrl("app1"), au.findUrlForWds(apps));
  }

  @Test
  void preferWdsAppOverCromwell() throws Exception {
    List<ListAppResponse> apps =
        new java.util.ArrayList<>(List.of(separatedWdsAppWdsAppName, separatedWorkflowsApp));

    permuteAndTest(apps, anticipatedWdsUrl("wds"));
  }

  @Test
  void throwIfBestAppNotReady() {
    List<ListAppResponse> apps =
        new java.util.ArrayList<>(
            List.of(otherNamedCromwellAppOlder, otherNamedCromwellAppProvisioning));
    // Shuffle to make sure the initial ordering isn't relevant:
    Collections.shuffle(apps);

    assertThrows(DependencyNotAvailableException.class, () -> au.findUrlForWds(apps));
  }

  @Test
  void throwIfBestAppHasNoWDS() {
    List<ListAppResponse> apps = new java.util.ArrayList<>(List.of(separatedWorkflowsApp));
    // Shuffle to make sure the initial ordering isn't relevant:
    Collections.shuffle(apps);

    assertThrows(DependencyNotAvailableException.class, () -> au.findUrlForWds(apps));
  }

  @Test
  void findCromwellRunnerUrlInCombinedResponse() throws Exception {
    List<ListAppResponse> apps = List.of(cromwellListAppResponse, cromwellRunnerAppResponse);

    AppUtils au = new AppUtils(leonardoServerConfiguration, imputationConfiguration);
    assertEquals(anticipatedCromwellUrl("cromwell-runner"), au.findUrlForCromwell(apps));
  }

  @Test
  void findCromwellUrlInCombinedWDSApp() throws Exception {
    List<ListAppResponse> apps = List.of(cromwellListAppResponse, separatedWdsAppWdsAppName);

    AppUtils au = new AppUtils(leonardoServerConfiguration, imputationConfiguration);
    assertEquals(anticipatedCromwellUrl("terra-app"), au.findUrlForCromwell(apps));
  }

  @Test
  void findCromwellRunnerUrlInCombinedWDSApp() throws Exception {
    List<ListAppResponse> apps = List.of(cromwellRunnerAppResponse, separatedWdsAppWdsAppName);

    AppUtils au = new AppUtils(leonardoServerConfiguration, imputationConfiguration);
    assertEquals(anticipatedCromwellUrl("cromwell-runner"), au.findUrlForCromwell(apps));
  }

  private void permuteAndTest(List<ListAppResponse> apps, String expectedUrl) throws Exception {
    int permutation = 0;
    do {
      Collections.rotate(apps, 1);
      assertEquals(expectedUrl, au.findUrlForWds(apps));
    } while (permutation++ < apps.size());
  }

  private String anticipatedWdsUrl(String appName) {
    return StringSubstitutor.replace(
        "https://lzblahblahblah.servicebus.windows.net/${appName}-${workspaceId}/wds",
        Map.of("workspaceId", workspaceId, "appName", appName));
  }

  private String anticipatedCromwellUrl(String appName) {
    return StringSubstitutor.replace(
        "https://lzblahblahblah.servicebus.windows.net/${appName}-${workspaceId}/cromwell",
        Map.of("workspaceId", workspaceId, "appName", appName));
  }

  public AppUtilsTest() throws IOException {
    // Use GSON instead of objectMapper because we want to simulate the JSON that comes back from
    // Leonardo.
    org.broadinstitute.dsde.workbench.client.leonardo.JSON.setGson(new Gson());
    combinedWdsInCromwellApp =
        ListAppResponse.fromJson(
            StringSubstitutor.replace(
                """
                {
                    "workspaceId": "${workspaceId}",
                    "cloudContext": {
                        "cloudProvider": "AZURE",
                        "cloudResource": "blah-blah-blah"
                    },
                    "kubernetesRuntimeConfig": {
                        "numNodes": 1,
                        "machineType": "Standard_A2_v2",
                        "autoscalingEnabled": false
                    },
                    "errors": [],
                    "status": "RUNNING",
                    "proxyUrls": {
                        "cbas": "https://lzblahblahblah.servicebus.windows.net/wds-${workspaceId}/cbas",
                        "cbas-ui": "https://lzblahblahblah.servicebus.windows.net/wds-${workspaceId}/",
                        "cromwell": "https://lzblahblahblah.servicebus.windows.net/wds-${workspaceId}/cromwell",
                        "wds": "https://lzblahblahblah.servicebus.windows.net/wds-${workspaceId}/wds"
                    },
                    "appName": "wds-${workspaceId}",
                    "appType": "CROMWELL",
                    "diskName": null,
                    "auditInfo": {
                        "creator": "me@broadinstitute.org",
                        "createdDate": "2023-02-09T16:01:36.660590Z",
                        "destroyedDate": null,
                        "dateAccessed": "2023-02-09T16:01:36.660590Z"
                    },
                    "accessScope": null,
                    "labels": {}
                }""",
                Map.of("workspaceId", workspaceId)));

    separatedWdsAppWdsAppName =
        ListAppResponse.fromJson(
            StringSubstitutor.replace(
                """
                    {
                        "workspaceId": "${workspaceId}",
                        "cloudContext": {
                            "cloudProvider": "AZURE",
                            "cloudResource": "blah-blah-blah"
                        },
                        "kubernetesRuntimeConfig": {
                            "numNodes": 1,
                            "machineType": "Standard_A2_v2",
                            "autoscalingEnabled": false
                        },
                        "errors": [],
                        "status": "RUNNING",
                        "proxyUrls": {
                            "wds": "https://lzblahblahblah.servicebus.windows.net/wds-${workspaceId}/wds"
                        },
                        "appName": "wds-${workspaceId}",
                        "appType": "WDS",
                        "diskName": null,
                        "auditInfo": {
                            "creator": "me@broadinstitute.org",
                            "createdDate": "2022-10-10T16:01:36.660590Z",
                            "destroyedDate": null,
                            "dateAccessed": "2023-02-09T16:01:36.660590Z"
                        },
                        "accessScope": null,
                        "labels": {}
                    }""",
                Map.of("workspaceId", workspaceId)));

    separatedWdsAppNotWdsAppName =
        ListAppResponse.fromJson(
            StringSubstitutor.replace(
                """
                                {
                                    "workspaceId": "${workspaceId}",
                                    "cloudContext": {
                                        "cloudProvider": "AZURE",
                                        "cloudResource": "blah-blah-blah"
                                    },
                                    "kubernetesRuntimeConfig": {
                                        "numNodes": 1,
                                        "machineType": "Standard_A2_v2",
                                        "autoscalingEnabled": false
                                    },
                                    "errors": [],
                                    "status": "RUNNING",
                                    "proxyUrls": {
                                        "wds": "https://lzblahblahblah.shouldneverwork.servicebus.windows.net/wds-${workspaceId}/wds"
                                    },
                                    "appName": "notwds-${workspaceId}",
                                    "appType": "WDS",
                                    "diskName": null,
                                    "auditInfo": {
                                        "creator": "me@broadinstitute.org",
                                        "createdDate": "2022-10-10T16:01:36.660590Z",
                                        "destroyedDate": null,
                                        "dateAccessed": "2023-02-09T16:01:36.660590Z"
                                    },
                                    "accessScope": null,
                                    "labels": {}
                                }""",
                Map.of("workspaceId", workspaceId)));

    cromwellListAppResponse =
        ListAppResponse.fromJson(
            StringSubstitutor.replace(
                """
                {
                    "workspaceId": "${workspaceId}",
                    "cloudContext": {
                        "cloudProvider": "AZURE",
                        "cloudResource": "blah-blah-blah"
                    },
                    "kubernetesRuntimeConfig": {
                        "numNodes": 1,
                        "machineType": "Standard_A2_v2",
                        "autoscalingEnabled": false
                    },
                    "errors": [],
                    "status": "RUNNING",
                    "proxyUrls": {
                        "cbas": "https://lzblahblahblah.servicebus.windows.net/terra-app-${workspaceId}/cbas",
                        "cbas-ui": "https://lzblahblahblah.servicebus.windows.net/terra-app-${workspaceId}/",
                        "cromwell": "https://lzblahblahblah.servicebus.windows.net/terra-app-${workspaceId}/cromwell"
                    },
                    "appName": "terra-app-${workspaceId}",
                    "appType": "CROMWELL",
                    "diskName": null,
                    "auditInfo": {
                        "creator": "me@broadinstitute.org",
                        "createdDate": "2023-02-09T16:01:36.660590Z",
                        "destroyedDate": null,
                        "dateAccessed": "2023-02-09T16:01:36.660590Z"
                    },
                    "accessScope": null,
                    "labels": {}
                }""",
                Map.of("workspaceId", workspaceId)));

    otherNamedCromwellAppOlder =
        ListAppResponse.fromJson(
            StringSubstitutor.replace(
                """
                {
                    "workspaceId": "${workspaceId}",
                    "cloudContext": {
                        "cloudProvider": "AZURE",
                        "cloudResource": "blah-blah-blah"
                    },
                    "kubernetesRuntimeConfig": {
                        "numNodes": 1,
                        "machineType": "Standard_A2_v2",
                        "autoscalingEnabled": false
                    },
                    "errors": [],
                    "status": "RUNNING",
                    "proxyUrls": {
                        "cbas": "https://lzblahblahblah.servicebus.windows.net/app2-${workspaceId}/cbas",
                        "cbas-ui": "https://lzblahblahblah.servicebus.windows.net/app2-${workspaceId}/",
                        "cromwell": "https://lzblahblahblah.servicebus.windows.net/app2-${workspaceId}/cromwell",
                        "wds": "https://lzblahblahblah.servicebus.windows.net/app2-${workspaceId}/wds"
                    },
                    "appName": "app2-${workspaceId}",
                    "appType": "CROMWELL",
                    "diskName": null,
                    "auditInfo": {
                        "creator": "me@broadinstitute.org",
                        "createdDate": "2022-10-10T16:01:36.660590Z",
                        "destroyedDate": null,
                        "dateAccessed": "2023-02-09T16:01:36.660590Z"
                    },
                    "accessScope": null,
                    "labels": {}
                }""",
                Map.of("workspaceId", workspaceId)));

    galaxyApp =
        ListAppResponse.fromJson(
            StringSubstitutor.replace(
                """
                {
                    "workspaceId": "${workspaceId}",
                    "cloudContext": {
                        "cloudProvider": "AZURE",
                        "cloudResource": "blah-blah-blah"
                    },
                    "kubernetesRuntimeConfig": {
                        "numNodes": 1,
                        "machineType": "Standard_A2_v2",
                        "autoscalingEnabled": false
                    },
                    "errors": [],
                    "status": "RUNNING",
                    "proxyUrls": {
                        "blah": "blah blah"
                    },
                    "appName": "galaxy-${workspaceId}",
                    "appType": "GALAXY",
                    "diskName": null,
                    "auditInfo": {
                        "creator": "me@broadinstitute.org",
                        "createdDate": "2023-02-09T16:01:36.660590Z",
                        "destroyedDate": null,
                        "dateAccessed": "2023-02-09T16:01:36.660590Z"
                    },
                    "accessScope": null,
                    "labels": {}
                }""",
                Map.of("workspaceId", workspaceId)));

    otherNamedCromwellApp =
        ListAppResponse.fromJson(
            StringSubstitutor.replace(
                """
                {
                    "workspaceId": "${workspaceId}",
                    "cloudContext": {
                        "cloudProvider": "AZURE",
                        "cloudResource": "blah-blah-blah"
                    },
                    "kubernetesRuntimeConfig": {
                        "numNodes": 1,
                        "machineType": "Standard_A2_v2",
                        "autoscalingEnabled": false
                    },
                    "errors": [],
                    "status": "RUNNING",
                    "proxyUrls": {
                        "cbas": "https://lzblahblahblah.servicebus.windows.net/app1-${workspaceId}/cbas",
                        "cbas-ui": "https://lzblahblahblah.servicebus.windows.net/app1-${workspaceId}/",
                        "cromwell": "https://lzblahblahblah.servicebus.windows.net/app1-${workspaceId}/cromwell",
                        "wds": "https://lzblahblahblah.servicebus.windows.net/app1-${workspaceId}/wds"
                    },
                    "appName": "app1-${workspaceId}",
                    "appType": "CROMWELL",
                    "diskName": null,
                    "auditInfo": {
                        "creator": "me@broadinstitute.org",
                        "createdDate": "2023-02-09T16:01:36.660590Z",
                        "destroyedDate": null,
                        "dateAccessed": "2023-02-09T16:01:36.660590Z"
                    },
                    "accessScope": null,
                    "labels": {}
                }""",
                Map.of("workspaceId", workspaceId)));

    otherNamedCromwellAppProvisioning =
        ListAppResponse.fromJson(
            StringSubstitutor.replace(
                """
                {
                    "workspaceId": "${workspaceId}",
                    "cloudContext": {
                        "cloudProvider": "AZURE",
                        "cloudResource": "blah-blah-blah"
                    },
                    "kubernetesRuntimeConfig": {
                        "numNodes": 1,
                        "machineType": "Standard_A2_v2",
                        "autoscalingEnabled": false
                    },
                    "errors": [],
                    "status": "PROVISIONING",
                    "proxyUrls": {
                        "cbas": "https://lzblahblahblah.servicebus.windows.net/app1-${workspaceId}/cbas",
                        "cbas-ui": "https://lzblahblahblah.servicebus.windows.net/app1-${workspaceId}/",
                        "cromwell": "https://lzblahblahblah.servicebus.windows.net/app1-${workspaceId}/cromwell",
                        "wds": "https://lzblahblahblah.servicebus.windows.net/app1-${workspaceId}/wds"
                    },
                    "appName": "app1-${workspaceId}",
                    "appType": "CROMWELL",
                    "diskName": null,
                    "auditInfo": {
                        "creator": "me@broadinstitute.org",
                        "createdDate": "2023-02-09T16:01:36.660590Z",
                        "destroyedDate": null,
                        "dateAccessed": "2023-02-09T16:01:36.660590Z"
                    },
                    "accessScope": null,
                    "labels": {}
                }""",
                Map.of("workspaceId", workspaceId)));

    separatedWorkflowsApp =
        ListAppResponse.fromJson(
            StringSubstitutor.replace(
                """
                {
                    "workspaceId": "${workspaceId}",
                    "cloudContext": {
                        "cloudProvider": "AZURE",
                        "cloudResource": "blah-blah-blah"
                    },
                    "kubernetesRuntimeConfig": {
                        "numNodes": 1,
                        "machineType": "Standard_A2_v2",
                        "autoscalingEnabled": false
                    },
                    "errors": [],
                    "status": "RUNNING",
                    "proxyUrls": {
                        "cbas": "https://lzblahblahblah.servicebus.windows.net/workflows-${workspaceId}/cbas",
                        "cbas-ui": "https://lzblahblahblah.servicebus.windows.net/workflows-${workspaceId}/",
                        "cromwell": "https://lzblahblahblah.servicebus.windows.net/workflows-${workspaceId}/cromwell"
                    },
                    "appName": "workflows-${workspaceId}",
                    "appType": "CROMWELL",
                    "diskName": null,
                    "auditInfo": {
                        "creator": "me@broadinstitute.org",
                        "createdDate": "2023-02-09T16:01:36.660590Z",
                        "destroyedDate": null,
                        "dateAccessed": "2023-02-09T16:01:36.660590Z"
                    },
                    "accessScope": null,
                    "labels": {}
                }""",
                Map.of("workspaceId", workspaceId)));

    cromwellRunnerAppResponse =
        ListAppResponse.fromJson(
            StringSubstitutor.replace(
                """
                {
                    "workspaceId": "${workspaceId}",
                    "cloudContext": {
                        "cloudProvider": "AZURE",
                        "cloudResource": "blah-blah-blah"
                    },
                    "kubernetesRuntimeConfig": {
                        "numNodes": 1,
                        "machineType": "Standard_A2_v2",
                        "autoscalingEnabled": false
                    },
                    "errors": [],
                    "status": "RUNNING",
                    "proxyUrls": {
                        "cromwell-runner": "https://lzblahblahblah.servicebus.windows.net/cromwell-runner-${workspaceId}/cromwell"
                    },
                    "appName": "cromwell-runner-${workspaceId}",
                    "appType": "CROMWELL_RUNNER_APP",
                    "diskName": null,
                    "auditInfo": {
                        "creator": "me@broadinstitute.org",
                        "createdDate": "2023-02-09T16:01:36.660590Z",
                        "destroyedDate": null,
                        "dateAccessed": "2023-02-09T16:01:36.660590Z"
                    },
                    "accessScope": null,
                    "labels": {}
                }""",
                Map.of("workspaceId", workspaceId)));
  }
}
