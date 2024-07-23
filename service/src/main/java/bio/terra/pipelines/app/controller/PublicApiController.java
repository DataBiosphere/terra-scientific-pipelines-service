package bio.terra.pipelines.app.controller;

import bio.terra.pipelines.app.configuration.internal.OidcConfiguration;
import bio.terra.pipelines.generated.api.PublicApi;
import bio.terra.pipelines.generated.model.ApiVersionProperties;
import bio.terra.pipelines.service.StatusService;
import jakarta.servlet.http.HttpServletResponse;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.util.UriComponentsBuilder;

@Controller
public class PublicApiController implements PublicApi {
  private final ApiVersionProperties versionProperties;
  private final StatusService statusService;
  private final OidcConfiguration oidcConfiguration;

  @Autowired
  public PublicApiController(
      ApiVersionProperties versionProperties,
      StatusService statusService,
      OidcConfiguration oidcConfiguration) {
    this.versionProperties = versionProperties;
    this.statusService = statusService;
    this.oidcConfiguration = oidcConfiguration;
    new ApiVersionProperties()
        .build(versionProperties.getBuild())
        .gitHash(versionProperties.getGitHash())
        .github(versionProperties.getGithub())
        .gitTag(versionProperties.getGitTag());
  }

  @Override
  public ResponseEntity<Void> getStatus() {
    return new ResponseEntity<>(
        statusService.getCurrentStatus() ? HttpStatus.OK : HttpStatus.SERVICE_UNAVAILABLE);
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

  //  @GetMapping(value = "/")
  //  public String index() {
  //    return "redirect:swagger-ui.html";
  //  }

  @GetMapping(value = "/swagger-ui.html")
  public String getSwagger(Model model) {
    model.addAttribute("clientId", oidcConfiguration.clientId());
    return "index";
  }

  @GetMapping(value = "/openapi.yml")
  public String getOpenApiYaml(Model model, HttpServletResponse response) {
    model.addAttribute("authorityEndpoint", getOidcMetadataEndpoint());
    // set CORS headers for the openapi.yml file to allow the central swagger-ui to fetch it
    response.setHeader("Access-Control-Allow-Origin", "*");
    response.setHeader("Access-Control-Allow-Methods", "GET, OPTIONS");
    response.setHeader("Access-Control-Allow-Headers", "Content-Type, api_key, Authorization");
    return "openapi";
  }

  private String getOidcMetadataEndpoint() {
    // parse the authority endpoint as a uri then append the oidc metadata path
    return UriComponentsBuilder.fromUriString(oidcConfiguration.authorityEndpoint())
        .path("/.well-known/openid-configuration")
        .toUriString();
  }
}
