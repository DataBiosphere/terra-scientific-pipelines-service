package bio.terra.pipelines.app.controller;

import bio.terra.pipelines.generated.api.PipelinesApi;
import bio.terra.pipelines.generated.model.ApiTspsPipeline;
import bio.terra.pipelines.generated.model.ApiTspsPipelinesGetResult;
import bio.terra.pipelines.service.pipelines.PipelinesService;
import bio.terra.pipelines.service.pipelines.model.Pipeline;
import java.util.List;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;

/** Main TSPS controller */
@Controller
public class PipelinesApiController implements PipelinesApi {
  private final PipelinesService pipelinesService;

  @Autowired
  public PipelinesApiController(PipelinesService pipelinesService) {
    this.pipelinesService = pipelinesService;
  }

  // -- Pipelines --

  @Override
  public ResponseEntity<ApiTspsPipelinesGetResult> getPipelines() {
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    ApiTspsPipelinesGetResult result = pipelinesToApi(pipelineList);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  static ApiTspsPipelinesGetResult pipelinesToApi(List<Pipeline> pipelineList) {
    ApiTspsPipelinesGetResult apiResult = new ApiTspsPipelinesGetResult();

    for (Pipeline pipeline : pipelineList) {
      var apiPipeline =
          new ApiTspsPipeline()
              .pipelineId(pipeline.getPipelineId())
              .displayName(pipeline.getDisplayName())
              .description(pipeline.getDescription());
      // TODO: is there a better function to use here? e.g. addPipelineItem()
      apiResult.add(0, apiPipeline);
    }

    return apiResult;
  }
}

// -- Scientific Pipelines Attribute Objects --
