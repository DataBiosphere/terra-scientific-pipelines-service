package bio.terra.pipelines.app.controller;

import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.generated.api.QuotasApi;
import bio.terra.pipelines.generated.model.ApiQuotaWithDetails;
import bio.terra.pipelines.service.QuotasService;
import io.swagger.annotations.Api;
import jakarta.servlet.http.HttpServletRequest;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;

/** Quotas controller */
@Controller
@Api(tags = {"quotas"})
public class QuotasController implements QuotasApi {
  private final SamConfiguration samConfiguration;
  private final SamUserFactory samUserFactory;
  private final HttpServletRequest request;
  private final QuotasService quotasService;

  private final Logger logger = LoggerFactory.getLogger(QuotasController.class);

  private SamUser getAuthenticatedInfo() {
    return samUserFactory.from(request, samConfiguration.baseUri());
  }

  QuotasController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      QuotasService quotasService) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.quotasService = quotasService;
  }

  /** Returns the quota usage for the authenticated user for the specified pipeline. */
  @Override
  public ResponseEntity<ApiQuotaWithDetails> getQuotaForPipeline(String pipelineName) {
    SamUser user = getAuthenticatedInfo();
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            user.getSubjectId(), validatedPipelineName);
    QuotaUnitsEnum quotaUnits = quotasService.getQuotaUnitsForPipeline(validatedPipelineName);

    return new ResponseEntity<>(quotasToApiQuotaWithDetails(userQuota, quotaUnits), HttpStatus.OK);
  }

  static ApiQuotaWithDetails quotasToApiQuotaWithDetails(
      UserQuota userQuota, QuotaUnitsEnum quotaUnits) {

    return new ApiQuotaWithDetails()
        .pipelineName(userQuota.getPipelineName().getValue())
        .quotaConsumed(userQuota.getQuotaConsumed())
        .quotaLimit(userQuota.getQuota())
        .quotaUnits(quotaUnits.getValue());
  }
}
