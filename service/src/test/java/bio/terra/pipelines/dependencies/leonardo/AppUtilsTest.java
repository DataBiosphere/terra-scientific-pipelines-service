package bio.terra.pipelines.dependencies.leonardo;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.app.configuration.external.LeonardoServerConfiguration;
import com.google.gson.Gson;
import java.io.IOException;
import java.util.*;
import org.apache.commons.text.StringSubstitutor;
import org.broadinstitute.dsde.workbench.client.leonardo.model.AppStatus;
import org.broadinstitute.dsde.workbench.client.leonardo.model.AppType;
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

  private final ListAppResponse separatedWdsAppWdsAppName;
  // does not have "wds-" app name and does not have anticipated wds url
  private final ListAppResponse separatedWdsAppNotWdsAppName;
  private final ListAppResponse cromwellListAppResponse;
  private final ListAppResponse cromwellListAppResponseProvisioning;
  private final ListAppResponse cromwellListAppResponseNoProxyUrlDeleted;
  private final ListAppResponse combinedWdsInCromwellApp;
  private final ListAppResponse otherNamedCromwellAppOlder;
  private final ListAppResponse galaxyApp;
  private final ListAppResponse otherNamedCromwellApp;
  private final ListAppResponse otherNamedCromwellAppProvisioning;
  private final ListAppResponse separatedWorkflowsApp;

  private final AppUtils au = new AppUtils(leonardoServerConfiguration);

  @Test
  void findWdsUrlInCombinedApp() {
    List<ListAppResponse> apps = List.of(combinedWdsInCromwellApp);

    AppUtils au = new AppUtils(leonardoServerConfiguration);
    assertEquals(anticipatedWdsUrl("wds"), au.findUrlForWds(apps, workspaceId));
  }

  @Test
  void findCromwellUrlInTerraApp() {
    List<ListAppResponse> apps = List.of(cromwellListAppResponse);

    AppUtils au = new AppUtils(leonardoServerConfiguration);
    assertEquals(anticipatedCbasUrl("terra-app"), au.findUrlForCbas(apps, workspaceId));
  }

  @Test
  void findWdsUrlInWdsApp() {
    List<ListAppResponse> apps = List.of(separatedWdsAppWdsAppName);

    AppUtils au = new AppUtils(leonardoServerConfiguration);
    assertEquals(anticipatedWdsUrl("wds"), au.findUrlForWds(apps, workspaceId));
  }

  @Test
  void findWdsAndCromwellInCombinedAppResponse() {
    List<ListAppResponse> apps = List.of(separatedWdsAppWdsAppName, cromwellListAppResponse);

    AppUtils au = new AppUtils(leonardoServerConfiguration);
    assertEquals(anticipatedWdsUrl("wds"), au.findUrlForWds(apps, workspaceId));
    assertEquals(anticipatedCbasUrl("terra-app"), au.findUrlForCbas(apps, workspaceId));
  }

  @Test
  void preferSpecificallyNamedApp() {
    List<ListAppResponse> apps =
        new java.util.ArrayList<>(List.of(separatedWdsAppWdsAppName, separatedWdsAppNotWdsAppName));
    permuteAndTest(apps, anticipatedWdsUrl("wds"));
  }

  @Test
  void notChooseGalaxyApp() {
    List<ListAppResponse> apps =
        new java.util.ArrayList<>(List.of(otherNamedCromwellApp, galaxyApp));
    // Shuffle to make sure the initial ordering isn't relevant:
    Collections.shuffle(apps);
    assertEquals(anticipatedWdsUrl("app1"), au.findUrlForWds(apps, workspaceId));
  }

  @Test
  void preferNewerCreatedApp() {
    List<ListAppResponse> apps =
        new java.util.ArrayList<>(List.of(otherNamedCromwellApp, otherNamedCromwellAppOlder));
    // Shuffle to make sure the initial ordering isn't relevant:
    Collections.shuffle(apps);
    assertEquals(anticipatedWdsUrl("app1"), au.findUrlForWds(apps, workspaceId));
  }

  @Test
  void preferWdsAppOverCromwell() {
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

    assertThrows(DependencyNotAvailableException.class, () -> au.findUrlForWds(apps, workspaceId));
  }

  @Test
  void throwIfBestAppHasNoWDS() {
    List<ListAppResponse> apps = new java.util.ArrayList<>(List.of(separatedWorkflowsApp));
    // Shuffle to make sure the initial ordering isn't relevant:
    Collections.shuffle(apps);

    assertThrows(DependencyNotAvailableException.class, () -> au.findUrlForWds(apps, workspaceId));
  }

  @Test
  void throwIfBestAppHasNoCromwell() {
    List<ListAppResponse> apps =
        new java.util.ArrayList<>(List.of(cromwellListAppResponseNoProxyUrlDeleted));
    // Shuffle to make sure the initial ordering isn't relevant:
    Collections.shuffle(apps);

    assertThrows(DependencyNotAvailableException.class, () -> au.findUrlForCbas(apps, workspaceId));
  }

  @Test
  void throwIfBestAppHasNoRunningCromwell() {
    List<ListAppResponse> apps =
        new java.util.ArrayList<>(List.of(cromwellListAppResponseProvisioning));
    // Shuffle to make sure the initial ordering isn't relevant:
    Collections.shuffle(apps);

    assertThrows(DependencyNotAvailableException.class, () -> au.findUrlForCbas(apps, workspaceId));
  }

  @Test
  void findCromwellUrlInCombinedWDSApp() {
    List<ListAppResponse> apps = List.of(cromwellListAppResponse, separatedWdsAppWdsAppName);

    AppUtils au = new AppUtils(leonardoServerConfiguration);
    assertEquals(anticipatedCbasUrl("terra-app"), au.findUrlForCbas(apps, workspaceId));
  }

  @Test
  void testFilterSuitableApps() {
    Set<AppStatus> healthyStates = EnumSet.of(AppStatus.RUNNING, AppStatus.STOPPED);
    List<AppType> appTypeList = List.of(AppType.CROMWELL, AppType.CROMWELL_RUNNER_APP);
    AppType appTypeToFilterOn = AppType.CROMWELL;

    // test an app that should pass the filter
    ListAppResponse goodApp = getPassingListAppResponse();
    assertTrue(
        au.passesSuitableFilter(
            goodApp, healthyStates, appTypeToFilterOn, appTypeList, workspaceId));

    // test an app that should fail the filter due to bad workspace id
    ListAppResponse badAppWorkspaceId =
        getPassingListAppResponse().workspaceId(UUID.randomUUID().toString());
    assertFalse(
        au.passesSuitableFilter(
            badAppWorkspaceId, healthyStates, appTypeToFilterOn, appTypeList, workspaceId));

    // test an app that should fail the filter due to app type
    ListAppResponse badAppAppType = getPassingListAppResponse().appType(AppType.WDS);
    assertFalse(
        au.passesSuitableFilter(
            badAppAppType, healthyStates, appTypeToFilterOn, appTypeList, workspaceId));

    // test an app that should fail the filter due to status
    ListAppResponse badAppStatus = getPassingListAppResponse().status(AppStatus.DELETED);
    assertFalse(
        au.passesSuitableFilter(
            badAppStatus, healthyStates, appTypeToFilterOn, appTypeList, workspaceId));
  }

  private ListAppResponse getPassingListAppResponse() {
    return new ListAppResponse()
        .workspaceId(workspaceId)
        .appType(AppType.CROMWELL)
        .status(AppStatus.RUNNING);
  }

  private void permuteAndTest(List<ListAppResponse> apps, String expectedUrl) {
    int permutation = 0;
    do {
      Collections.rotate(apps, 1);
      assertEquals(expectedUrl, au.findUrlForWds(apps, workspaceId));
    } while (permutation++ < apps.size());
  }

  private String anticipatedWdsUrl(String appName) {
    return StringSubstitutor.replace(
        "https://lzblahblahblah.servicebus.windows.net/${appName}-${workspaceId}/wds",
        Map.of("workspaceId", workspaceId, "appName", appName));
  }

  private String anticipatedCbasUrl(String appName) {
    return StringSubstitutor.replace(
        "https://lzblahblahblah.servicebus.windows.net/${appName}-${workspaceId}/cbas",
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

    cromwellListAppResponseNoProxyUrlDeleted =
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
                                "status": "DELETED",
                                "proxyUrls": {
                                    "cbas": "https://lzblahblahblah.servicebus.windows.net/terra-app-${workspaceId}/cbas",
                                    "cbas-ui": "https://lzblahblahblah.servicebus.windows.net/terra-app-${workspaceId}/"
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

    cromwellListAppResponseProvisioning =
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
  }
}
