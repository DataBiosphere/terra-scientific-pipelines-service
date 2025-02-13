package bio.terra.pipelines.app.configuration.external;

import org.springframework.boot.context.properties.ConfigurationProperties;

@ConfigurationProperties(prefix = "teaspoons.ingress")
public class IngressConfiguration {

  /** Fully-qualified domain name. The base URL this instance can be accessed at. */
  private String domainName;

  public String getDomainName() {
    return domainName;
  }

  public void setDomainName(String domainName) {
    this.domainName = domainName;
  }
}
