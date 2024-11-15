package bio.terra.pipelines.testutils;

import bio.terra.common.iam.BearerToken;
import bio.terra.common.iam.SamUser;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.UserQuota;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.ObjectWriter;
import com.fasterxml.jackson.databind.SerializationFeature;
import java.util.UUID;
import org.springframework.test.web.servlet.request.MockHttpServletRequestBuilder;

/**
 * A collection of utilities and constants useful for MockMVC-based tests. This style of tests lets
 * us test controller-layer code (request/response parsing, authz, and validation) without actually
 * spinning up a local server.
 */
public class MockMvcUtils {
  public static final String AUTH_HEADER = "Authorization";

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
    // needed because we currently allow amorphous objects currently as pipelineInputs
    mapper.configure(SerializationFeature.FAIL_ON_EMPTY_BEANS, false);
    ObjectWriter ow = mapper.writer();
    return ow.writeValueAsString(obj);
  }

  // Common test constants
  public static final SamUser TEST_SAM_USER =
      new SamUser(
          "test@email",
          UUID.randomUUID().toString(),
          new BearerToken(UUID.randomUUID().toString()));

  public static final UUID TEST_WORKSPACE_UUID =
      UUID.fromString("94fd136b-1234-1234-1234-76d8a2811066");
  public static final String TEST_WORKSPACE_BILLING_PROJECT = "testTerraProject";
  public static final String TEST_WORKSPACE_NAME = "testTerraWorkspaceName";
  public static final String TEST_WORKSPACE_STORAGE_CONTAINER_NAME = "test-bucket-name";
  public static final String TEST_WORKSPACE_GOOGLE_PROJECT = "testGoogleProject";
  public static final String TEST_WDL_METHOD_VERSION = "0.12.1";
  // using this function to build a pipeline with a value set for the id field.  Normally this would
  // be populated
  // by calling `save()` from the repository but since these tests mock that out, we have to set the
  // value of id
  // ourselves.
  public static Pipeline getTestPipeline() {
    Pipeline testPipeline =
        new Pipeline(
            PipelinesEnum.ARRAY_IMPUTATION,
            0,
            "displayName",
            "description",
            "pipelineType",
            "wdlUrl",
            "wdlMethodName",
            TEST_WDL_METHOD_VERSION,
            TEST_WORKSPACE_UUID,
            TEST_WORKSPACE_BILLING_PROJECT,
            TEST_WORKSPACE_NAME,
            TEST_WORKSPACE_STORAGE_CONTAINER_NAME,
            TEST_WORKSPACE_GOOGLE_PROJECT,
            TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST,
            TestUtils.TEST_PIPELINE_OUTPUTS_DEFINITION_LIST);
    testPipeline.setId(2L);
    return testPipeline;
  }

  public static final UserQuota TEST_USER_QUOTA_1 =
      new UserQuota(PipelinesEnum.ARRAY_IMPUTATION, TEST_SAM_USER.getSubjectId(), 1000, 10);
}
