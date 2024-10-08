package bio.terra.pipelines.common.components;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.app.components.TemplateResolvers;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.thymeleaf.templatemode.TemplateMode;
import org.thymeleaf.templateresolver.ClassLoaderTemplateResolver;

class TemplateResolversTest extends BaseEmbeddedDbTest {

  @Autowired TemplateResolvers templateResolvers;

  @Test
  void secondaryTemplateResolver() {
    ClassLoaderTemplateResolver secondaryTemplateResolver =
        templateResolvers.secondaryTemplateResolver();
    assertNotNull(secondaryTemplateResolver);
    assertEquals("static/", secondaryTemplateResolver.getPrefix());
    assertEquals(".yml", secondaryTemplateResolver.getSuffix());
    assertEquals(TemplateMode.TEXT, secondaryTemplateResolver.getTemplateMode());
    assertEquals("UTF-8", secondaryTemplateResolver.getCharacterEncoding());
    assertEquals(1, secondaryTemplateResolver.getOrder());
    assertTrue(secondaryTemplateResolver.getCheckExistence());
  }
}
