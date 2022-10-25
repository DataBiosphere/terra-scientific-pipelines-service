package bio.terra.pipelines.app.controller;

import bio.terra.pipelines.generated.api.TspsApi;
import bio.terra.pipelines.service.pao.PipelinesService;
import bio.terra.pipelines.service.pao.model.Pipeline;
import java.util.List;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;

/** Main TSPS controller */
@Controller
public class TspsApiController implements TspsApi {
  private final PipelinesService pipelinesService;

  @Autowired
  public TspsApiController(PipelinesService pipelinesService) {
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
