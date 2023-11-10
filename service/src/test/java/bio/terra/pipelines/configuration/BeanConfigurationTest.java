package bio.terra.pipelines.configuration;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import bio.terra.common.exception.UnauthorizedException;
import bio.terra.pipelines.app.configuration.internal.BeanConfiguration;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.http.HttpHeaders;
import org.springframework.mock.web.MockHttpServletRequest;

class BeanConfigurationTest {
  private final BeanConfiguration beanConfiguration = new BeanConfiguration(null);
  private MockHttpServletRequest request;

  @BeforeEach
  void init() {
    request = new MockHttpServletRequest();
  }

  @Test
  void bearerTokenSuccess() {
    String token = "token";
    request.addHeader(HttpHeaders.AUTHORIZATION, "Bearer " + token);
    assertEquals(token, beanConfiguration.bearerToken(request).getToken());
  }

  @Test
  void noBearerToken() {
    assertThrows(UnauthorizedException.class, () -> beanConfiguration.bearerToken(request));
  }
}
