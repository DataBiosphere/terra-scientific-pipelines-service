package bio.terra.pipelines.app.controller;

import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineQuota;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.generated.api.PipelinesApi;
import bio.terra.pipelines.generated.model.*;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.service.QuotasService;
import io.swagger.annotations.Api;
import jakarta.servlet.http.HttpServletRequest;
import java.util.List;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.PathVariable;

/** Pipelines controller */
@Controller
@Api(tags = {"pipelines"})
public class PipelinesApiController implements PipelinesApi {
  private final SamConfiguration samConfiguration;
  private final SamUserFactory samUserFactory;
  private final HttpServletRequest request;
  private final PipelinesService pipelinesService;
  private final SamService samService;
  private final QuotasService quotasService;

  @Autowired
  public PipelinesApiController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      PipelinesService pipelinesService,
      SamService samService,
      QuotasService quotasService) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.pipelinesService = pipelinesService;
    this.samService = samService;
    this.quotasService = quotasService;
  }

  private static final Logger logger = LoggerFactory.getLogger(PipelinesApiController.class);

  private SamUser getAuthenticatedInfo() {
    return samUserFactory.from(request, samConfiguration.baseUri());
  }

  // -- Pipelines --

  @Override
  public ResponseEntity<ApiGetPipelinesResult> getPipelines() {
    SamUser authedUser = getAuthenticatedInfo();
    boolean showHiddenPipelines =
        samService.isAdmin(authedUser.getEmail(), authedUser.getBearerToken().getToken());
    List<Pipeline> pipelineList = pipelinesService.getPipelines(showHiddenPipelines);
    ApiGetPipelinesResult result = pipelinesToApi(pipelineList);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiPipelineWithDetails> getPipelineDetails(
      @PathVariable("pipelineName") String pipelineName, ApiGetPipelineDetailsRequestBody body) {
    SamUser authedUser = getAuthenticatedInfo();
    boolean showHiddenPipelines =
        samService.isAdmin(authedUser.getEmail(), authedUser.getBearerToken().getToken());
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);

    // Get the pipeline version from the request body, if it exists
    Integer pipelineVersion = body == null ? null : body.getPipelineVersion();
    Pipeline pipelineInfo =
        pipelinesService.getPipeline(validatedPipelineName, pipelineVersion, showHiddenPipelines);

    // Fetch the quota settings to attach to the pipeline details
    PipelineQuota pipelineQuota = quotasService.getPipelineQuota(validatedPipelineName);
    ApiPipelineQuota apiPipelineQuota = pipelineQuotaToApi(pipelineQuota);

    ApiPipelineWithDetails result = pipelineWithDetailsToApi(pipelineInfo);
    result.setPipelineQuota(apiPipelineQuota);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  static ApiGetPipelinesResult pipelinesToApi(List<Pipeline> pipelineList) {
    ApiGetPipelinesResult apiResult = new ApiGetPipelinesResult();

    apiResult.setResults(pipelineList.stream().map(PipelinesApiController::pipelineToApi).toList());

    return apiResult;
  }

  static ApiPipelineQuota pipelineQuotaToApi(PipelineQuota pipelineQuota) {
    return new ApiPipelineQuota()
        .pipelineName(pipelineQuota.getPipelineName().getValue())
        .defaultQuota(pipelineQuota.getDefaultQuota())
        .minQuotaConsumed(pipelineQuota.getMinQuotaConsumed())
        .quotaUnits(pipelineQuota.getQuotaUnits().toString());
  }

  static ApiPipelineWithDetails pipelineWithDetailsToApi(Pipeline pipelineInfo) {
    ApiPipelineUserProvidedInputDefinitions inputs = new ApiPipelineUserProvidedInputDefinitions();
    inputs.addAll(
        pipelineInfo.getPipelineInputDefinitions().stream()
            .filter(PipelineInputDefinition::isUserProvided)
            .map(
                input ->
                    new ApiPipelineUserProvidedInputDefinition()
                        .name(input.getName())
                        .displayName(input.getDisplayName())
                        .description(input.getDescription())
                        .type(input.getType().toString())
                        .isRequired(input.isRequired())
                        .defaultValue(input.getDefaultValue())
                        .minValue(input.getMinValue())
                        .maxValue(input.getMaxValue())
                        .fileSuffix(input.getFileSuffix()))
            .toList());
    ApiPipelineOutputDefinitions outputs = new ApiPipelineOutputDefinitions();
    outputs.addAll(
        pipelineInfo.getPipelineOutputDefinitions().stream()
            .map(
                output ->
                    new ApiPipelineOutputDefinition()
                        .name(output.getName())
                        .displayName(output.getDisplayName())
                        .description(output.getDescription())
                        .type(output.getType().toString()))
            .toList());
    return new ApiPipelineWithDetails()
        .pipelineName(pipelineInfo.getName().getValue())
        .displayName(pipelineInfo.getDisplayName())
        .pipelineVersion(pipelineInfo.getVersion())
        .description(pipelineInfo.getDescription())
        .type(pipelineInfo.getPipelineType())
        .inputs(inputs)
        .outputs(outputs);
  }

  static ApiPipeline pipelineToApi(Pipeline pipelineInfo) {
    return new ApiPipeline()
        .pipelineName(pipelineInfo.getName().getValue())
        .displayName(pipelineInfo.getDisplayName())
        .pipelineVersion(pipelineInfo.getVersion())
        .description(pipelineInfo.getDescription());
  }
}
