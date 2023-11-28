package bio.terra.pipelines.service;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoServiceException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import java.util.Collections;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

/** Service to encapsulate logic used to run an imputation pipeline */
@Service
public class ImputationService {
  private static final Logger logger = LoggerFactory.getLogger(ImputationService.class);
  private LeonardoService leonardoService;
  private SamService samService;
  private WdsService wdsService;
  private ImputationConfiguration imputationConfiguration;

  @Autowired
  ImputationService(
      LeonardoService leonardoService,
      SamService samService,
      WdsService wdsService,
      ImputationConfiguration imputationConfiguration) {
    this.leonardoService = leonardoService;
    this.samService = samService;
    this.wdsService = wdsService;
    this.imputationConfiguration = imputationConfiguration;
  }

  public List<ListAppResponse> queryForWorkspaceApps() {
    String workspaceId = imputationConfiguration.workspaceId();
    try {
      List<ListAppResponse> getAppsResponse =
          leonardoService.getApps(workspaceId, samService.getTspsServiceAccountToken(), false);

      logger.info(
          "GetAppsResponse for workspace id {}: {}",
          imputationConfiguration.workspaceId(),
          getAppsResponse);

      String wdsUri =
          leonardoService.getWdsUrlFromApps(
              workspaceId, samService.getTspsServiceAccountToken(), false);
      logger.info("Wds uri for workspace id {}: {}", imputationConfiguration.workspaceId(), wdsUri);
      logger.info(
          "Wds health: {}",
          wdsService.checkHealth(wdsUri, samService.getTspsServiceAccountToken()));

      logger.info(
          "Wds schema: {}",
          wdsService.querySchema(
              wdsUri,
              samService.getTspsServiceAccountToken(),
              imputationConfiguration.workspaceId()));

      return getAppsResponse;
    } catch (LeonardoServiceException e) {
      logger.error("Get Apps called for workspace id {} failed", workspaceId);
      return Collections.emptyList();
    } catch (WdsServiceException e) {
      logger.error("Calls to Wds for workspace id {} failed", workspaceId);
      return Collections.emptyList();
    }
  }
}
