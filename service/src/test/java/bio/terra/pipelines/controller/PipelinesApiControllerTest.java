package bio.terra.pipelines.controller;

import static bio.terra.pipelines.testutils.MockMvcUtils.getTestPipeline;
import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.*;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.post;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.common.iam.BearerTokenFactory;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.app.controller.GlobalExceptionHandler;
import bio.terra.pipelines.app.controller.PipelinesApiController;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineQuota;
import bio.terra.pipelines.db.exception.InvalidPipelineException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.generated.model.*;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.testutils.MockMvcUtils;
import bio.terra.pipelines.testutils.TestUtils;
import com.fasterxml.jackson.databind.ObjectMapper;
import jakarta.servlet.http.HttpServletRequest;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.http.MediaType;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.bean.override.mockito.MockitoBean;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.MvcResult;

@ContextConfiguration(classes = {PipelinesApiController.class, GlobalExceptionHandler.class})
@WebMvcTest
class PipelinesApiControllerTest {
  @MockitoBean PipelinesService pipelinesServiceMock;
  @MockitoBean QuotasService quotasServiceMock;
  @MockitoBean SamUserFactory samUserFactoryMock;
  @MockitoBean BearerTokenFactory bearerTokenFactory;
  @MockitoBean SamConfiguration samConfiguration;
  @MockitoBean SamService samService;

  @Autowired private MockMvc mockMvc;

  private final List<Pipeline> testPipelineList =
      List.of(TestUtils.TEST_PIPELINE_1, TestUtils.TEST_PIPELINE_2);
  private final SamUser testUser = MockMvcUtils.TEST_SAM_USER;

  @BeforeEach
  void beforeEach() {
    when(samConfiguration.baseUri()).thenReturn("baseSamUri");
    when(samUserFactoryMock.from(any(HttpServletRequest.class), eq("baseSamUri")))
        .thenReturn(testUser);
    when(pipelinesServiceMock.getPipeline(any(PipelinesEnum.class), anyInt()))
        .thenReturn(getTestPipeline());
    when(quotasServiceMock.getPipelineQuota(any(PipelinesEnum.class)))
        .thenReturn(TestUtils.TEST_PIPELINE_QUOTA_1);
  }

  // getPipeline tests

  @Test
  void getPipelinesOk() throws Exception {
    when(pipelinesServiceMock.getPipelines()).thenReturn(testPipelineList);

    MvcResult result =
        mockMvc
            .perform(get("/api/pipelines/v1"))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiGetPipelinesResult response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiGetPipelinesResult.class);

    assertEquals(testPipelineList.size(), response.getResults().size());
  }

  @Test
  void getPipelineDetailsOkNoVersion() throws Exception {
    PipelinesEnum pipelineNameEnum = TestUtils.TEST_PIPELINE_1.getName();
    String pipelineName = pipelineNameEnum.getValue();

    when(pipelinesServiceMock.getPipeline(pipelineNameEnum, null))
        .thenReturn(TestUtils.TEST_PIPELINE_1);

    MvcResult result =
        mockMvc
            .perform(
                post("/api/pipelines/v1/" + pipelineName)
                    .contentType(MediaType.APPLICATION_JSON)
                    .content("{}"))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiPipelineWithDetails response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiPipelineWithDetails.class);

    assertEquals(pipelineName, response.getPipelineName());
    assertEquals(TestUtils.TEST_PIPELINE_VERSION_1, response.getPipelineVersion());
    assertEquals(TestUtils.TEST_PIPELINE_1.getDescription(), response.getDescription());
    assertEquals(TestUtils.TEST_PIPELINE_1.getDisplayName(), response.getDisplayName());
    assertEquals(TestUtils.TEST_PIPELINE_1.getPipelineType(), response.getType());
    assertEquals(TestUtils.TEST_PIPELINE_1.getVersion(), response.getPipelineVersion());

    // check that the response only includes user-provided inputs
    assertEquals(
        TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST.stream()
            .filter(PipelineInputDefinition::isUserProvided)
            .toList()
            .size(),
        response.getInputs().size());
    for (ApiPipelineUserProvidedInputDefinition p : response.getInputs()) {
      // find the matching input definition in test pipeline inputs list and check if it's user
      // provided
      assertTrue(
          TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST.stream()
              .anyMatch(i -> i.getName().equals(p.getName()) && i.isUserProvided()));
    }
  }

  @Test
  void getPipelineDetailsOkWithVersion() throws Exception {
    PipelinesEnum pipelineNameEnum = TestUtils.TEST_PIPELINE_1.getName();
    String pipelineName = pipelineNameEnum.getValue();

    when(pipelinesServiceMock.getPipeline(pipelineNameEnum, 3))
        .thenReturn(TestUtils.TEST_PIPELINE_1);

    MvcResult result =
        mockMvc
            .perform(
                post("/api/pipelines/v1/" + pipelineName)
                    .contentType(MediaType.APPLICATION_JSON)
                    .content("{\"pipelineVersion\":\"%s\"}".formatted(3)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiPipelineWithDetails response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiPipelineWithDetails.class);

    assertEquals(pipelineName, response.getPipelineName());
    assertEquals(TestUtils.TEST_PIPELINE_1.getDescription(), response.getDescription());
    assertEquals(TestUtils.TEST_PIPELINE_1.getDisplayName(), response.getDisplayName());
    assertEquals(TestUtils.TEST_PIPELINE_1.getPipelineType(), response.getType());
    assertEquals(TestUtils.TEST_PIPELINE_1.getVersion(), response.getPipelineVersion());

    // check that the response only includes user-provided inputs
    assertEquals(
        TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST.stream()
            .filter(PipelineInputDefinition::isUserProvided)
            .toList()
            .size(),
        response.getInputs().size());
    for (ApiPipelineUserProvidedInputDefinition p : response.getInputs()) {
      // find the matching input definition in test pipeline inputs list and check if it's user
      // provided
      assertTrue(
          TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST.stream()
              .anyMatch(i -> i.getName().equals(p.getName()) && i.isUserProvided()));
    }
  }

  @Test
  void getPipelineDetailsIncludesQuota() throws Exception {
    PipelinesEnum pipelineNameEnum = TestUtils.TEST_PIPELINE_1.getName();
    String pipelineName = pipelineNameEnum.getValue();

    // Mocks
    when(pipelinesServiceMock.getPipeline(pipelineNameEnum, null))
        .thenReturn(TestUtils.TEST_PIPELINE_1);
    PipelineQuota testQuota = new PipelineQuota(pipelineNameEnum, 100, 15, QuotaUnitsEnum.SAMPLES);
    when(quotasServiceMock.getPipelineQuota(pipelineNameEnum)).thenReturn(testQuota);

    MvcResult result =
        mockMvc
            .perform(
                post("/api/pipelines/v1/" + pipelineName)
                    .contentType(MediaType.APPLICATION_JSON)
                    .content("{}"))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiPipelineWithDetails response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiPipelineWithDetails.class);

    assertNotNull(response.getPipelineQuota());
    assertEquals(pipelineName, response.getPipelineQuota().getPipelineName());
    assertEquals(100, response.getPipelineQuota().getDefaultQuota());
    assertEquals(15, response.getPipelineQuota().getMinQuotaConsumed());
    assertEquals("SAMPLES", response.getPipelineQuota().getQuotaUnits());
  }

  @Test
  void getPipelineCaseInsensitive() throws Exception {
    String pipelineName = "aRrAy_ImpuTatioN";
    PipelinesEnum pipelineNameEnum = PipelinesEnum.ARRAY_IMPUTATION;

    when(pipelinesServiceMock.getPipeline(pipelineNameEnum, null))
        .thenReturn(TestUtils.TEST_PIPELINE_1);

    MvcResult result =
        mockMvc
            .perform(
                post("/api/pipelines/v1/" + pipelineName)
                    .contentType(MediaType.APPLICATION_JSON)
                    .content("{}"))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiPipelineWithDetails response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiPipelineWithDetails.class);

    assertEquals(pipelineNameEnum.getValue(), response.getPipelineName());
  }

  @Test
  void getPipelineBadPipeline() throws Exception {
    String pipelineName = "bad-pipeline-id";

    mockMvc
        .perform(
            post("/api/pipelines/v1/" + pipelineName)
                .contentType(MediaType.APPLICATION_JSON)
                .content("{}"))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(InvalidPipelineException.class, result.getResolvedException()));
  }

  @Test
  void getPipelines500ErrorGenericSupportResponse() throws Exception {
    when(pipelinesServiceMock.getPipelines()).thenThrow(new RuntimeException("test exception"));

    MvcResult result =
        mockMvc
            .perform(get("/api/pipelines/v1"))
            .andExpect(status().is5xxServerError())
            .andReturn();

    ApiErrorReport response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiErrorReport.class);

    assertEquals(
        "Internal server error. Please contact support if this problem persists.",
        response.getMessage());
  }
}
