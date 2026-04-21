package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.app.configuration.internal.OidcConfiguration;
import bio.terra.pipelines.generated.api.PublicApi;
import bio.terra.pipelines.generated.model.ApiVersionProperties;
import bio.terra.pipelines.service.StatusService;
import jakarta.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.io.DefaultResourceLoader;
import org.springframework.core.io.Resource;
import org.springframework.core.io.ResourceLoader;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;
import org.springframework.ui.Model;
import org.springframework.util.FileCopyUtils;
import org.springframework.web.bind.annotation.GetMapping;

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
  public ResponseEntity<String> getDocs(String docType) {
    String resourcePath;
    if (docType.equals("termsOfService")) {
      resourcePath = "classpath:documents/TermsOfService.md";
    } else if (docType.equals("acceptableUsePolicy")) {
      resourcePath = "classpath:documents/AcceptableUse.md";
    } else {
      throw new BadRequestException("%s is not a valid docType".formatted(docType));
    }
    try {
      ResourceLoader resourceLoader = new DefaultResourceLoader();
      Resource resource = resourceLoader.getResource(resourcePath);
      String response =
          new String(
              FileCopyUtils.copyToByteArray(resource.getInputStream()), StandardCharsets.UTF_8);
      return new ResponseEntity<>(response, HttpStatus.OK);
    } catch (IOException e) {
      throw new InternalServerErrorException("Failed to load document: " + docType, e);
    }
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

  private static final String CSP_HEADER_NAME = "Content-Security-Policy";
  private static final String CSP_HEADER_CONTENTS =
      "script-src 'self' 'unsafe-inline'; img-src 'self' data:; style-src 'self' 'unsafe-inline'; form-action 'none';";

  @GetMapping({"/", "/index.html", "/swagger-ui.html"})
  public String getSwagger(Model model, HttpServletResponse response) {
    response.setHeader(CSP_HEADER_NAME, CSP_HEADER_CONTENTS);

    model.addAttribute("clientId", oidcConfiguration.clientId());
    return "index";
  }

  @GetMapping(value = "/openapi.yml")
  public String getOpenApiYaml(Model model, HttpServletResponse response) {
    response.setHeader(CSP_HEADER_NAME, CSP_HEADER_CONTENTS);

    model.addAttribute("authorityEndpoint", oidcConfiguration.authorityEndpoint());
    model.addAttribute("tokenEndpoint", oidcConfiguration.tokenEndpoint());
    return "openapi";
  }
}
