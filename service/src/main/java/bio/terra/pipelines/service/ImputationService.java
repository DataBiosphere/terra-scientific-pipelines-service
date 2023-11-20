package bio.terra.pipelines.service;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoServiceException;
import bio.terra.pipelines.dependencies.sam.SamService;
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
  private ImputationConfiguration imputationConfiguration;

  @Autowired
  ImputationService(
      LeonardoService leonardoService,
      SamService samService,
      ImputationConfiguration imputationConfiguration) {
    this.leonardoService = leonardoService;
    this.samService = samService;
    this.imputationConfiguration = imputationConfiguration;
  }

  public List<ListAppResponse> queryForWorkspaceApps() {
    try {
      List<ListAppResponse> getAppsResponse =
          leonardoService.getApps(
              imputationConfiguration.workspaceId(),
              samService.getTspsServiceAccountToken(),
              false);

      logger.info(
          "GetAppsResponse for workspace id {}: {}",
          imputationConfiguration.workspaceId(),
          getAppsResponse.toString());
      return getAppsResponse;
    } catch (LeonardoServiceException e) {
      logger.error(
          "Get Apps called for workspace id {} failed", imputationConfiguration.workspaceId());
      return null;
    }
  }
}
