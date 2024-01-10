package bio.terra.pipelines.controller;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.post;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.common.exception.MissingRequiredFieldException;
import bio.terra.common.iam.BearerTokenFactory;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.app.controller.PipelinesApiController;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.exception.InvalidPipelineException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.generated.model.ApiCreateJobRequestBody;
import bio.terra.pipelines.generated.model.ApiGetPipelinesResult;
import bio.terra.pipelines.generated.model.ApiPipeline;
import bio.terra.pipelines.service.ImputationService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.testutils.MockMvcUtils;
import bio.terra.pipelines.testutils.TestUtils;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.util.List;
import java.util.UUID;
import javax.servlet.http.HttpServletRequest;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.http.MediaType;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.MvcResult;

@ContextConfiguration(classes = PipelinesApiController.class)
@WebMvcTest
class PipelinesApiControllerTest {
  @MockBean PipelinesService pipelinesServiceMock;
  @MockBean StairwayJobService stairwayJobServiceMock;
  @MockBean SamUserFactory samUserFactoryMock;
  @MockBean BearerTokenFactory bearerTokenFactory;
  @MockBean SamConfiguration samConfiguration;
  @MockBean SamService samService;
  @MockBean ImputationService imputationService;

  @Autowired private MockMvc mockMvc;

  private final List<Pipeline> testPipelineList =
      List.of(TestUtils.TEST_PIPELINE_1, TestUtils.TEST_PIPELINE_2);
  private final SamUser testUser = MockMvcUtils.TEST_SAM_USER;
  private final String testPipelineVersion = TestUtils.TEST_PIPELINE_VERSION_1;
  private final Object testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;

  @BeforeEach
  void beforeEach() {
    when(samUserFactoryMock.from(any(HttpServletRequest.class), any())).thenReturn(testUser);
    when(imputationService.queryForWorkspaceApps()).thenReturn(null);
  }

  @Test
  void testGetPipelinesOk() throws Exception {
    when(pipelinesServiceMock.getPipelines()).thenReturn(testPipelineList);

    MvcResult result =
        mockMvc
            .perform(get("/api/pipelines/v1alpha1"))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiGetPipelinesResult response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiGetPipelinesResult.class);

    assertEquals(testPipelineList.size(), response.size());
  }

  @Test
  void getPipelineOk() throws Exception {
    String pipelineId = TestUtils.TEST_PIPELINE_1.getPipelineId();
    when(pipelinesServiceMock.getPipeline(pipelineId)).thenReturn(TestUtils.TEST_PIPELINE_1);

    MvcResult result =
        mockMvc
            .perform(get("/api/pipelines/v1alpha1/" + pipelineId))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiPipeline response =
        new ObjectMapper().readValue(result.getResponse().getContentAsString(), ApiPipeline.class);

    assertEquals(pipelineId, response.getPipelineId());
  }

  @Test
  void testCreateJobImputationPipeline() throws Exception {
    String pipelineId = PipelinesEnum.IMPUTATION.getValue();
    // This makes the body of the post... which is a lot for very little
    ApiCreateJobRequestBody postBody =
        new ApiCreateJobRequestBody()
            .pipelineVersion(testPipelineVersion)
            .pipelineInputs(testPipelineInputs);
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    UUID jobId = UUID.randomUUID(); // newJobId

    // the mocks
    when(imputationService.createImputationJob(
            testUser.getSubjectId(), testPipelineVersion, testPipelineInputs))
        .thenReturn(jobId);

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelines/v1alpha1/%s", pipelineId))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    assertTrue(result.getResponse().getContentAsString().contains(jobId.toString()));
  }

  @Test
  void testCreateJobBadPipeline() throws Exception {
    String pipelineId = "bad-pipeline-id";

    // This makes the body of the post... which is a lot for very little
    ApiCreateJobRequestBody postBody =
        new ApiCreateJobRequestBody()
            .pipelineVersion(testPipelineVersion)
            .pipelineInputs(testPipelineInputs);
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // no mocks since this should throw on validatePipelineId()
    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1alpha1/%s", pipelineId))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isNotFound())
        .andExpect(
            result ->
                assertTrue(result.getResolvedException() instanceof InvalidPipelineException));
  }

  @Test
  void testCreateImputationStairwayError() throws Exception {
    String pipelineId = PipelinesEnum.IMPUTATION.getValue();

    // This makes the body of the post... which is a lot for very little
    ApiCreateJobRequestBody postBody =
        new ApiCreateJobRequestBody()
            .pipelineVersion(testPipelineVersion)
            .pipelineInputs(testPipelineInputs);
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // the mocks - one error that can happen is a MissingRequiredFieldException from Stairway
    when(imputationService.createImputationJob(
            testUser.getSubjectId(), testPipelineVersion, testPipelineInputs))
        .thenThrow(new MissingRequiredFieldException("some message"));

    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1alpha1/%s", pipelineId))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isInternalServerError())
        .andExpect(
            result ->
                assertTrue(result.getResolvedException() instanceof MissingRequiredFieldException));
  }
}
