package bio.terra.pipelines.testutils;

import bio.terra.pipelines.service.model.Pipeline;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.ObjectWriter;
import com.fasterxml.jackson.databind.SerializationFeature;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.test.web.servlet.request.MockHttpServletRequestBuilder;

/**
 * A collection of utilities and constants useful for MockMVC-based tests. This style of tests lets
 * us test controller-layer code (request/response parsing, authz, and validation) without actually
 * spinning up a local server.
 */
public class MockMvcUtils {
  public static final String AUTH_HEADER = "Authorization";
  private static final Logger logger = LoggerFactory.getLogger(MockMvcUtils.class);

  public static MockHttpServletRequestBuilder addAuth(MockHttpServletRequestBuilder request) {
    return request.header(AUTH_HEADER, "Bearer ThisIsNotARealBearerToken");
  }

  public static MockHttpServletRequestBuilder addJsonContentType(
      MockHttpServletRequestBuilder request) {
    return request.contentType("application/json");
  }

  public static String convertToJsonString(Object obj) throws JsonProcessingException {
    ObjectMapper mapper = new ObjectMapper();
    mapper.configure(SerializationFeature.WRAP_ROOT_VALUE, false);
    ObjectWriter ow = mapper.writer();
    return ow.writeValueAsString(obj);
  }

  // Pipelines test constants
  public static final String TEST_PIPELINE_ID_1 = "test-pipeline-id-1";
  public static final String TEST_PIPELINE_NAME_1 = "Test Pipeline Name One";
  public static final String TEST_PIPELINE_DESCRIPTION_1 = "Test Pipeline Description One";
  public static final String TEST_PIPELINE_ID_2 = "test-pipeline-id-2";
  public static final String TEST_PIPELINE_NAME_2 = "Test Pipeline Name Two";
  public static final String TEST_PIPELINE_DESCRIPTION_2 = "Test Pipeline Description Two";
  public static final Pipeline TEST_PIPELINE_1 =
      new Pipeline(TEST_PIPELINE_ID_1, TEST_PIPELINE_NAME_1, TEST_PIPELINE_DESCRIPTION_1);
  public static final Pipeline TEST_PIPELINE_2 =
      new Pipeline(TEST_PIPELINE_ID_2, TEST_PIPELINE_NAME_2, TEST_PIPELINE_DESCRIPTION_2);

  public static final String TEST_USER_ID_1 = "test-user-id-1";
  public static final String TEST_USER_ID_2 = "test-user-id-2";
}
