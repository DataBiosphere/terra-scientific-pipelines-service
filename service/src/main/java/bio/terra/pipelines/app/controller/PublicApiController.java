package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.BadRequestException;
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
    String response;
    if (docType.equals("termsOfService")) {
      response =
"""
# Data Science Services at Broad Clinical Laboratories - Product Specific Terms of Service

## 1. Introduction

Thank you for using Data Science Services, an application provided by Broad Clinical Laboratories, LLC ("BCL" or "we"). Data Science Services, or the application it runs on, the Terra Scientific Pipelines Service, offers registered users of the application ("Users") an automated and scalable platform for Users to perform genomics analyses in a secure environment.

These Product Terms of Service (the "Product Terms") supplement the Broad Clinical Laboratories Terms and Conditions of Service set forth above (the "Terms") for your use of Data Science Services as operated by BCL. Any defined terms not defined herein carry the definition set forth in the Terms.

## 2. Your Use of Data Science Services

- **Own Risk:** Data Science Services are provided as-is and at your own risk.
- **Account and Security:** Users must create a Terra account and accept the Terra Terms of Service (or execute a separate license agreement to Terra). Data Science Services can be operated via a web portal User Interface (UI), with a command-line interface (CLI) tool called terralab, or directly via the service's APIs. You or your organization are responsible for the activity that happens on or through your account and for securing your log in information.
- **Monitoring and Improvement:** We may monitor your conduct, your use of Data Science Services, and your submitted data to ensure security, fix bugs, comply with the Terms, including the Product Terms, and Applicable Law, and improve the Service.
- **Publishing:** If publishing results that derived from data processed through Data Science Services, please cite the service as follows (whitepaper citation to be added):
  > Data Science Services at Broad Clinical Laboratories. (date of access: Year, Month Day). Pipeline name (version). https://services.terra.bio/
  >
  > e.g. Data Science Services at Broad Clinical Laboratories. (2026, Feb 18). All of Us + AnVIL Array Imputation (v1). https://services.terra.bio/

## 3. Submitted Data and Service Operation

- **Service Operation:** Data Science Services is a highly scalable, cloud-based service, maintained as NIST-800-53 Rev 5 Moderate and FedRAMP-compliant. The service runs genomic analyses, such as imputation, on the User's behalf against data that is protected and not directly accessible to the User.
- **Your Submitted Data:** You retain ownership of any intellectual property rights that you hold in that submitted data. Your submitted data will never be shared through Data Science Services with any other User without your permission. By submitting data to Data Science Services, you represent to us that you own, or have all necessary rights, to provide that data to Data Science Services.
- **Data Retention:** Once your imputed data is generated, you will have a limited time to download the data before the files are deleted.
- **Permitted Uses:** By submitting your data, you grant us a non-exclusive, royalty-free, worldwide license to use that submitted data for the limited purposes of: (a) processing your submitted data and returning the analysis result to you (i.e., the imputed data), (b) providing support, (c) using certain metadata about the service (but not your genomic data itself) to maintain, test, and improve Data Science Services, and (d) to comply with Applicable Law.

## 5. Pricing and Quota

Data Science Services utilizes a per-user quota system to manage costs and ensure fair access.

- **Quota System:** Initial access to pipelines (e.g., Public Beta, Friends & Family) may be offered free of charge but with a capped quota.
  - Creating alternate accounts in an attempt to gain more free quota is prohibited and may result in our suspending or permanently disabling your access to Data Science Services.
- **Payment:** For-profit or high-volume usage will require user payment and involve purchasing quota. See Section 4 of the Terms for additional information.

## 6. Modifying and Terminating Use

- **Changes to Data Science Services:** While we have no obligation to maintain or update Data Science Services, we are constantly changing and improving the service. We may make performance or security improvements, change functionalities or features, or make changes to comply with Applicable Law or to prevent illegal activities on, or abuse of, our systems. Your continued use of Data Science Services constitutes acceptance of any updates to the Terms, including these Product Terms.
- **Notification of Changes to the Service:** We publish any user-facing changes in our documentation. There may be times when we will need to make changes to Data Science Services without giving notice. We will seek to limit such instances to cases where we need to take action to ensure the security and operability of the service, prevent abuse or where we must act to comply with Applicable Laws.
""";
    } else if (docType.equals("acceptableUsePolicy")) {
      response =
"""
## Acceptable Uses

You may use Data Science Services for all lawful purposes, except as set forth below.

## Unacceptable Uses

You may not use Data Science Services for any of the following:

- **Sensitive Personal Data:** Submitting data that includes Protected Health Information (PHI) or special category personal data (as defined under GDPR, UK GDPR, or other applicable privacy laws), unless such data is fully deidentified or anonymized prior to submission or you first enter into a data processing agreement with BCL to provide such data in compliance with Applicable Laws. See Section 8 of the Terms for additional information.
- **Illegal or Harmful Activity:**
  - Any use in violation of Applicable Laws.
  - To try and gain unauthorized access or to disrupt any service, device, data, account, or network.
  - In a way that harms Data Science Services or impairs its use by other Users.
- **Data Security Breaches:**
  - Attempting to extract protected data from the service.
  - Processing national security information classified under US or other applicable law.
  - Browsing, searching, revealing, or retrieving protected data or information, or in any other way disclosing information, for someone who does not have authority to access that information.
  - Uploading any content that contains a software virus or malware, such as a Trojan horse or any other computer codes, files, or programs that may alter, damage, or interrupt the functioning of the service.
- **Circumvention:** Establishing any unauthorized interfaces between systems, networks, and applications owned or controlled by BCL or The Broad Institute and any other systems, network, and applications, whether owned and controlled by BCL, The Broad Institute, or a third party.
- **Assistance:** To assist, encourage, or induce any User or third party to do any of the above.
""";
    } else {
      throw new BadRequestException("%s is not a valid docType".formatted(docType));
    }
    return new ResponseEntity<>(response, HttpStatus.OK);
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
