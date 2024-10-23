package bio.terra.pipelines.app.controller;

import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.generated.api.PipelinesApi;
import bio.terra.pipelines.generated.model.ApiGetPipelinesResult;
import bio.terra.pipelines.generated.model.ApiPipeline;
import bio.terra.pipelines.generated.model.ApiPipelineUserProvidedInputDefinition;
import bio.terra.pipelines.generated.model.ApiPipelineUserProvidedInputDefinitions;
import bio.terra.pipelines.generated.model.ApiPipelineWithDetails;
import bio.terra.pipelines.service.PipelinesService;
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

  @Autowired
  public PipelinesApiController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      PipelinesService pipelinesService) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.pipelinesService = pipelinesService;
  }

  private static final Logger logger = LoggerFactory.getLogger(PipelinesApiController.class);

  private SamUser getAuthenticatedInfo() {
    return samUserFactory.from(request, samConfiguration.baseUri());
  }

  // -- Pipelines --

  @Override
  public ResponseEntity<ApiGetPipelinesResult> getPipelines() {
    getAuthenticatedInfo();
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    ApiGetPipelinesResult result = pipelinesToApi(pipelineList);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiPipelineWithDetails> getPipelineDetails(
      @PathVariable("pipelineName") String pipelineName) {
    getAuthenticatedInfo();
    PipelinesEnum validatedPipelineName =
        PipelineApiUtils.validatePipelineName(pipelineName, logger);
    Pipeline pipelineInfo = pipelinesService.getPipeline(validatedPipelineName);
    ApiPipelineWithDetails result = pipelineWithDetailsToApi(pipelineInfo);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  static ApiGetPipelinesResult pipelinesToApi(List<Pipeline> pipelineList) {
    ApiGetPipelinesResult apiResult = new ApiGetPipelinesResult();

    apiResult.setResults(pipelineList.stream().map(PipelinesApiController::pipelineToApi).toList());

    return apiResult;
  }

  static ApiPipelineWithDetails pipelineWithDetailsToApi(Pipeline pipelineInfo) {
    ApiPipelineUserProvidedInputDefinitions inputs = new ApiPipelineUserProvidedInputDefinitions();
    inputs.addAll(
        pipelineInfo.getPipelineInputDefinitions().stream()
            .filter(PipelineInputDefinition::getUserProvided)
            .map(
                input ->
                    new ApiPipelineUserProvidedInputDefinition()
                        .name(input.getName())
                        .type(input.getType().toString())
                        .isRequired(input.getIsRequired()))
            .toList());
    return new ApiPipelineWithDetails()
        .pipelineName(pipelineInfo.getName().getValue())
        .displayName(pipelineInfo.getDisplayName())
        .description(pipelineInfo.getDescription())
        .type(pipelineInfo.getPipelineType())
        .inputs(inputs);
  }

  static ApiPipeline pipelineToApi(Pipeline pipelineInfo) {
    return new ApiPipeline()
        .pipelineName(pipelineInfo.getName().getValue())
        .displayName(pipelineInfo.getDisplayName())
        .description(pipelineInfo.getDescription());
  }
}
