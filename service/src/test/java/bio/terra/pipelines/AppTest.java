package bio.terra.pipelines;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.thymeleaf.templatemode.TemplateMode;
import org.thymeleaf.templateresolver.ClassLoaderTemplateResolver;

class AppTest extends BaseEmbeddedDbTest {

  @Autowired App app;

  @Test
  void secondaryTemplateResolver() {
    ClassLoaderTemplateResolver secondaryTemplateResolver = app.secondaryTemplateResolver();
    assertNotNull(secondaryTemplateResolver);
    assertEquals("static/", secondaryTemplateResolver.getPrefix());
    assertEquals(".yml", secondaryTemplateResolver.getSuffix());
    assertEquals(TemplateMode.TEXT, secondaryTemplateResolver.getTemplateMode());
    assertEquals("UTF-8", secondaryTemplateResolver.getCharacterEncoding());
    assertEquals(1, secondaryTemplateResolver.getOrder());
    assertTrue(secondaryTemplateResolver.getCheckExistence());
  }
}
