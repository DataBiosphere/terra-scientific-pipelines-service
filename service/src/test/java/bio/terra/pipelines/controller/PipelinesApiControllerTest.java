package bio.terra.pipelines.controller;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.common.iam.BearerTokenFactory;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.app.controller.PipelinesApiController;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.generated.model.ApiGetPipelinesResult;
import bio.terra.pipelines.generated.model.ApiPipeline;
import bio.terra.pipelines.service.ImputationService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.testutils.MockMvcUtils;
import bio.terra.pipelines.testutils.TestUtils;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.util.List;
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

  @BeforeEach
  void beforeEach() {
    when(samUserFactoryMock.from(any(HttpServletRequest.class), any())).thenReturn(testUser);
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

  // TODO come back and fix these
  //  @Test
  //  void testCreateJobGoodPipeline() throws Exception {
  //    // This makes the body of the post... which is a lot for very little
  //    ApiCreateJobRequestBody postBody =
  //            new
  // ApiCreateJobRequestBody().pipelineVersion(pipelineVersion).pipelineInputs(pipelineInputs);
  //    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);
  //
  //    UUID fakeJobId = UUID.randomUUID();
  //
  //    // the mocks
  //    when(jobsServiceMock.createJob(
  //            testUser.getSubjectId(), pipelineId, pipelineVersion, pipelineInputs))
  //            .thenReturn(fakeJobId);
  //    when(pipelinesServiceMock.pipelineExists(pipelineId)).thenReturn(true);
  //
  //    // the crafting the expected response json
  //    Map<String, UUID> expectedResponseMap = new HashMap<>();
  //    expectedResponseMap.put("jobId", fakeJobId);
  //    String expectedResponseJson = MockMvcUtils.convertToJsonString(expectedResponseMap);
  //    mockMvc
  //            .perform(
  //                    post(String.format("/api/jobs/v1alpha1/%s", pipelineId))
  //                            .contentType(MediaType.APPLICATION_JSON)
  //                            .content(postBodyAsJson))
  //            .andExpect(status().isOk())
  //            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
  //            .andExpect(content().string(expectedResponseJson));
  //  }
  //
  //  @Test
  //  void testCreateJobBadPipeline() throws Exception {
  //    // This makes the body of the post... which is a lot for very little
  //    ApiPostJobRequestBody postBody =
  //            new
  // ApiPostJobRequestBody().pipelineVersion(pipelineVersion).pipelineInputs(pipelineInputs);
  //    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);
  //
  //    // the mocks
  //    when(pipelinesServiceMock.pipelineExists(pipelineId)).thenReturn(false);
  //
  //    mockMvc
  //            .perform(
  //                    post(String.format("/api/jobs/v1alpha1/%s", pipelineId))
  //                            .contentType(MediaType.APPLICATION_JSON)
  //                            .content(postBodyAsJson))
  //            .andExpect(status().isNotFound())
  //            .andExpect(
  //                    result ->
  //                            assertTrue(result.getResolvedException() instanceof
  // PipelineNotFoundException));
  //  }
  //
  //  @Test
  //  void testCreateJobWriteFail() throws Exception {
  //    // This makes the body of the post... which is a lot for very little
  //    ApiPostJobRequestBody postBody =
  //            new
  // ApiPostJobRequestBody().pipelineVersion(pipelineVersion).pipelineInputs(pipelineInputs);
  //    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);
  //
  //    // the mocks - if createJob repeatedly fails to write to the database, it returns null
  //    when(pipelinesServiceMock.pipelineExists(pipelineId)).thenReturn(true);
  //    when(jobsServiceMock.createJob(
  //            testUser.getSubjectId(), pipelineId, pipelineVersion, pipelineInputs))
  //            .thenReturn(null);
  //
  //    mockMvc
  //            .perform(
  //                    post(String.format("/api/jobs/v1alpha1/%s", pipelineId))
  //                            .contentType(MediaType.APPLICATION_JSON)
  //                            .content(postBodyAsJson))
  //            .andExpect(status().isInternalServerError())
  //            .andExpect(result -> assertTrue(result.getResolvedException() instanceof
  // ApiException));
  //  }
}
