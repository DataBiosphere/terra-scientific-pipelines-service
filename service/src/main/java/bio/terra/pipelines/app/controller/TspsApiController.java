package bio.terra.pipelines.app.controller;

import bio.terra.pipelines.generated.api.TspsApi;
import bio.terra.pipelines.service.pao.PaoService;
import java.util.UUID;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;

/** Main TSPS controller */
@Controller
public class TspsApiController implements TspsApi {
  private final PaoService paoService;

  @Autowired
  public TspsApiController(PaoService paoService) {
    this.paoService = paoService;
  }

  // -- Policy Queries --
  // TODO: PF-1733 Next step is to add group membership constraint

  // -- Policy Attribute Objects --
  @Override
  public ResponseEntity<Void> deletePao(UUID objectId) {
    paoService.deletePao(objectId);
    return new ResponseEntity<>(HttpStatus.NO_CONTENT);
  }
}

// -- Scientific Pipelines Attribute Objects --
