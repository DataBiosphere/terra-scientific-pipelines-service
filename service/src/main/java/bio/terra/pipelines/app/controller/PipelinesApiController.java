package bio.terra.pipelines.app.controller;

import bio.terra.pipelines.generated.api.PipelinesApi;
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
  public ResponseEntity<List> getPipelines() {
    List<Pipeline> pipelineList = pipelinesService.getPipelines();
    //    ApiPipelinesGetResult result = pipelineList;

    return new ResponseEntity<>(pipelineList, HttpStatus.OK);
  }
}

// -- Scientific Pipelines Attribute Objects --
