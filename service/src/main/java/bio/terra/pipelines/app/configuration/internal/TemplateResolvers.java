package bio.terra.pipelines.app.configuration.internal;

import org.springframework.context.annotation.Bean;
import org.springframework.stereotype.Component;
import org.thymeleaf.templatemode.TemplateMode;
import org.thymeleaf.templateresolver.ClassLoaderTemplateResolver;

@Component
public class TemplateResolvers {

  /**
   * This bean is used to resolve the location of the Thymeleaf templates that are used to generate
   * the OpenAPI documentation, i.e. static/openapi.yml as referenced in templates/index.html. On
   * the other hand, the default resolver is used to resolve the location of the Swagger UI
   * index.html file in templates/.
   */
  @Bean
  public ClassLoaderTemplateResolver secondaryTemplateResolver() {
    ClassLoaderTemplateResolver secondaryTemplateResolver = new ClassLoaderTemplateResolver();
    secondaryTemplateResolver.setPrefix("static/");
    secondaryTemplateResolver.setSuffix(".yml");
    secondaryTemplateResolver.setTemplateMode(TemplateMode.TEXT);
    secondaryTemplateResolver.setCharacterEncoding("UTF-8");
    secondaryTemplateResolver.setOrder(1);
    secondaryTemplateResolver.setCheckExistence(true);

    return secondaryTemplateResolver;
  }
}
