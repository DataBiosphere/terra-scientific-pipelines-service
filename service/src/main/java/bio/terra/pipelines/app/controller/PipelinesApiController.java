package bio.terra.pipelines.app.controller;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.generated.api.PipelinesApi;
import bio.terra.pipelines.generated.model.ApiPipeline;
import bio.terra.pipelines.generated.model.ApiPipelinesGetResult;
import bio.terra.pipelines.service.PipelinesService;
import io.swagger.annotations.Api;
import java.util.List;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.PathVariable;

/** Pipelines controller */
@Controller
@Api(tags = {"pipelines"})
public class PipelinesApiController implements PipelinesApi {
  private final PipelinesService pipelinesService;

  @Autowired
  public PipelinesApiController(PipelinesService pipelinesService) {
    this.pipelinesService = pipelinesService;
  }

  // -- Pipelines --

  @Override
  public ResponseEntity<ApiPipelinesGetResult> getPipelines() {
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    ApiPipelinesGetResult result = pipelinesToApi(pipelineList);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiPipeline> getPipeline(@PathVariable("pipelineId") String pipelineId) {
    Pipeline pipelineInfo = pipelinesService.getImputationPipelineViaFlight();
    ApiPipeline result = pipelineToApi(pipelineInfo);

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

  static ApiPipeline pipelineToApi(Pipeline pipelineInfo) {
    ApiPipeline apiPipeline =
        new ApiPipeline()
            .pipelineId(pipelineInfo.getPipelineId())
            .displayName(pipelineInfo.getDisplayName())
            .description(pipelineInfo.getDescription());

    return apiPipeline;
  }
}
