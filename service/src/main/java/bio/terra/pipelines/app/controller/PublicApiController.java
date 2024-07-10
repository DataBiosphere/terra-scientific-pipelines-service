package bio.terra.pipelines.app.controller;

import bio.terra.pipelines.generated.api.PublicApi;
import bio.terra.pipelines.generated.model.ApiVersionProperties;
import bio.terra.pipelines.service.StatusService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;

@Controller
public class PublicApiController implements PublicApi {
  private final ApiVersionProperties versionProperties;
  private final StatusService statusService;

  @Autowired
  public PublicApiController(ApiVersionProperties versionProperties, StatusService statusService) {
    this.versionProperties = versionProperties;
    this.statusService = statusService;
    new ApiVersionProperties()
        .build(versionProperties.getBuild())
        .gitHash(versionProperties.getGitHash())
        .github(versionProperties.getGithub())
        .gitTag(versionProperties.getGitTag());
  }

  @Override
  public ResponseEntity<Void> getStatus() {
    return new ResponseEntity<>(
        statusService.getCurrentStatus() ? HttpStatus.OK : HttpStatus.INTERNAL_SERVER_ERROR);
  }

  @Override
  public ResponseEntity<ApiVersionProperties> getVersion() {
    // these copy shenanigans are because versionProperties comes from spring config
    // and is actually a proxy and using the instance directly in a http response includes all the
    // proxy fields that no one wants to see
    return ResponseEntity.ok(
        new ApiVersionProperties()
            .build(versionProperties.getBuild())
            .gitHash(versionProperties.getGitHash())
            .github(versionProperties.getGithub())
            .gitTag(versionProperties.getGitTag()));
  }
}
