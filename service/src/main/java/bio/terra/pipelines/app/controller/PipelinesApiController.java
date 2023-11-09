package bio.terra.pipelines.app.controller;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.generated.api.PipelinesApi;
import bio.terra.pipelines.generated.model.ApiPipeline;
import bio.terra.pipelines.generated.model.ApiPipelinesGetResult;
import bio.terra.pipelines.service.PipelinesService;
import io.swagger.annotations.Api;
import java.util.List;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;

/** Pipelines controller */
@Controller
@Api(tags = {"pipelines"})
public class PipelinesApiController implements PipelinesApi {
  private final PipelinesService pipelinesService;
  private final LeonardoService leonardoService;

  @Autowired
  public PipelinesApiController(
      PipelinesService pipelinesService, LeonardoService leonardoService) {
    this.pipelinesService = pipelinesService;
    this.leonardoService = leonardoService;
  }

  private static final Logger logger = LoggerFactory.getLogger(PipelinesApiController.class);

  // -- Pipelines --

  @Override
  public ResponseEntity<ApiPipelinesGetResult> getPipelines() {
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    ApiPipelinesGetResult result = pipelinesToApi(pipelineList);
    // remove this before merge, was used to test we could actually get a result
    logger.info(leonardoService.getApps(false).toString());

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  static ApiPipelinesGetResult pipelinesToApi(List<Pipeline> pipelineList) {
    ApiPipelinesGetResult apiResult = new ApiPipelinesGetResult();

    for (Pipeline pipeline : pipelineList) {
      var apiPipeline =
          new ApiPipeline()
              .pipelineId(pipeline.getPipelineId())
              .displayName(pipeline.getDisplayName())
              .description(pipeline.getDescription());

      apiResult.add(apiPipeline);
    }

    return apiResult;
  }
}
