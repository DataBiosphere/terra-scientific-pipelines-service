package bio.terra.pipelines.dependencies.leonardo;

import bio.terra.pipelines.app.configuration.external.LeonardoServerConfiguration;
import java.time.OffsetDateTime;
import java.util.*;
import java.util.stream.Collectors;
import org.broadinstitute.dsde.workbench.client.leonardo.model.AppStatus;
import org.broadinstitute.dsde.workbench.client.leonardo.model.AppType;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

@Component
public class AppUtils {

  private final LeonardoServerConfiguration leonardoServerConfiguration;

  private static final Logger logger = LoggerFactory.getLogger(AppUtils.class);

  public AppUtils(LeonardoServerConfiguration leonardoServerConfiguration) {
    this.leonardoServerConfiguration = leonardoServerConfiguration;
  }

  int appComparisonFunction(
      ListAppResponse appToCompareA,
      ListAppResponse appToCompareB,
      List<AppType> appTypeList,
      String workspaceId) {
    // First criteria: Prefer apps with the expected app type.
    // NB: Negative because lower index is better
    int appTypeScoreA = -appTypeList.indexOf(appToCompareA.getAppType());
    int appTypeScoreB = -appTypeList.indexOf(appToCompareB.getAppType());
    if (appTypeScoreA != appTypeScoreB) {
      return appTypeScoreA - appTypeScoreB;
    }
    // If there is a WDS app type present, do this check; does not apply to cromwell app-types
    if (appToCompareA.getAppType() == AppType.WDS || appToCompareB.getAppType() == AppType.WDS) {
      // Second criteria: Prefer apps with the expected app type name
      int nameScoreA =
          Objects.equals(appToCompareA.getAppName(), "wds-%s".formatted(workspaceId)) ? 1 : 0;
      int nameScoreB =
          Objects.equals(appToCompareB.getAppName(), "wds-%s".formatted(workspaceId)) ? 1 : 0;
      if (nameScoreA != nameScoreB) {
        return nameScoreA - nameScoreB;
      }
    }

    // Third criteria: tie-break on whichever is older
    return OffsetDateTime.parse(appToCompareA.getAuditInfo().getCreatedDate())
        .compareTo(OffsetDateTime.parse(appToCompareB.getAuditInfo().getCreatedDate()));
  }

  /**
   * Invokes logic to determine the appropriate app for WDS and CROMWELL. If app is not running, a
   * URL will not be present, in this case we return empty string Note: This logic is similar to how
   * DataTable finds WDS app in Terra UI
   *
   * <p>(<a
   * href="https://github.com/DataBiosphere/terra-ui/blob/ac13bdf3954788ca7c8fd27b8fd4cfc755f150ff/src/libs/ajax/data-table-providers/WdsDataTableProvider.ts#L94-L147">...</a>)
   */
  ListAppResponse findBestAppForAppType(
      List<ListAppResponse> apps, AppType appType, String workspaceId)
      throws DependencyNotAvailableException {
    // WDS and Cromwell apps look for Kubernetes deployment statuses (such as RUNNING or
    // PROVISIONING), expressed by
    // Leo
    // See here for specific enumerations --
    // https://github.com/DataBiosphere/leonardo/blob/develop/core/src/main/scala/org/broadinstitute/dsde/workbench/leonardo/kubernetesModels.scala
    // look explicitly for a RUNNING app named 'wds-${app.workspaceId}' or
    // 'terra-app-<random_uuid>' -- if app is healthy and
    // running, there should only be one app RUNNING
    // an app may be in the 'PROVISIONING', 'STOPPED', 'STOPPING', which can still be deemed as an
    // OK state for Leonardo apps
    List<AppType> appTypeList;

    // WDS and CROMWELL apps will get their proxy urls from a different group of app types (located
    // in application.yml). Because of this, we need to identify the app type we need a url for, and
    // get the url from the correct group of relevant apps.
    if (appType.equals(AppType.WDS)) {
      appTypeList = leonardoServerConfiguration.wdsAppTypeNames();
    } else {
      appTypeList = leonardoServerConfiguration.cbasAppTypeNames();
    }

    Set<AppStatus> healthyStates =
        EnumSet.of(
            AppStatus.RUNNING,
            AppStatus.PROVISIONING,
            AppStatus.STARTING,
            AppStatus.STOPPED,
            AppStatus.STOPPING);

    List<ListAppResponse> suitableApps =
        apps.stream()
            .filter(
                app -> passesSuitableFilter(app, healthyStates, appType, appTypeList, workspaceId))
            .toList();

    return suitableApps.stream()
        .max(
            (appToCompareA, appToCompareB) ->
                this.appComparisonFunction(appToCompareA, appToCompareB, appTypeList, workspaceId))
        .orElseThrow(
            () ->
                new DependencyNotAvailableException(
                    "%s".formatted(appType.toString()),
                    "No suitable, healthy app found for %s (out of %s total apps in this workspace)"
                        .formatted(appType.toString(), apps.size())));
  }

  boolean passesSuitableFilter(
      ListAppResponse app,
      Set<AppStatus> healthyStates,
      AppType appType,
      List<AppType> appTypeList,
      String workspaceId) {
    var appMatchesWorkspaceId = Objects.equals(app.getWorkspaceId(), workspaceId);
    if (!appMatchesWorkspaceId) {
      logger.info(
          "Not using app {} for {} because it is in workspace {}, not {}",
          app.getAppName(),
          appType,
          app.getWorkspaceId(),
          workspaceId);
    }
    var appTypeListContainsApp = appTypeList.contains(app.getAppType());
    if (!appTypeListContainsApp) {
      logger.info(
          "Not using app {} for {} because it is of type {}, not one of {}",
          app.getAppName(),
          appType,
          app.getAppType(),
          appTypeList);
    }
    var isAppHealthy = healthyStates.contains(app.getStatus());
    if (!isAppHealthy) {
      logger.info(
          "Not using app {} for {} because it is in state {}, not one of {}",
          app.getAppName(),
          appType,
          app.getStatus(),
          healthyStates);
    }

    return appMatchesWorkspaceId && appTypeListContainsApp && isAppHealthy;
  }

  public String findUrlForWds(List<ListAppResponse> apps, String workspaceId)
      throws DependencyNotAvailableException {
    ListAppResponse foundApp = findBestAppForAppType(apps, AppType.WDS, workspaceId);
    @SuppressWarnings("unchecked")
    Map<String, String> proxyUrls = (foundApp.getProxyUrls());
    if (proxyUrls != null && foundApp.getStatus() == AppStatus.RUNNING) {
      return Optional.ofNullable(proxyUrls.get("wds"))
          .orElseThrow(
              () ->
                  new DependencyNotAvailableException(
                      "WDS",
                      "WDS proxy URL not found in %s app (available proxy URLs: %s)"
                          .formatted(
                              foundApp.getAppName(),
                              proxyUrls.keySet().stream()
                                  .sorted()
                                  .collect(Collectors.joining(", ")))));
    }

    throw new DependencyNotAvailableException(
        "WDS",
        "WDS in %s app not ready (%s)".formatted(foundApp.getAppName(), foundApp.getStatus()));
  }

  public String findUrlForCbas(List<ListAppResponse> apps, String workspaceId)
      throws DependencyNotAvailableException {
    ListAppResponse foundApp = findBestAppForAppType(apps, AppType.CROMWELL, workspaceId);

    String proxyUrl = "cbas";
    // find proper proxy for cromwell app type
    @SuppressWarnings("unchecked")
    Map<String, String> proxyUrls = foundApp.getProxyUrls();
    if (proxyUrls != null && foundApp.getStatus() == AppStatus.RUNNING) {
      return Optional.ofNullable(proxyUrls.get(proxyUrl))
          .orElseThrow(
              () ->
                  new DependencyNotAvailableException(
                      "Cromwell",
                      "CBAS proxy URL not found in %s app (available proxy URLs: %s)"
                          .formatted(
                              foundApp.getAppName(),
                              proxyUrls.keySet().stream()
                                  .sorted()
                                  .collect(Collectors.joining(", ")))));
    }

    throw new DependencyNotAvailableException(
        "Cromwell",
        "CBAS in %s app not ready (%s)".formatted(foundApp.getAppName(), foundApp.getStatus()));
  }
}
