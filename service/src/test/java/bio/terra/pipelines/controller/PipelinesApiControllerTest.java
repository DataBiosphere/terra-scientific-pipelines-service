package bio.terra.pipelines.controller;

import static org.junit.jupiter.api.Assertions.*;
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
import bio.terra.pipelines.app.controller.GlobalExceptionHandler;
import bio.terra.pipelines.app.controller.PipelinesApiController;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.exception.InvalidPipelineException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.dependencies.stairway.exception.InternalStairwayException;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJob;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.generated.model.*;
import bio.terra.pipelines.service.ImputationService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.testutils.MockMvcUtils;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import com.fasterxml.jackson.core.JsonProcessingException;
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

@ContextConfiguration(classes = {PipelinesApiController.class, GlobalExceptionHandler.class})
@WebMvcTest
class PipelinesApiControllerTest {
  @MockBean PipelinesService pipelinesServiceMock;
  @MockBean JobService jobServiceMock;
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
  private final UUID newJobId = TestUtils.TEST_NEW_UUID;

  @BeforeEach
  void beforeEach() {
    when(samUserFactoryMock.from(any(HttpServletRequest.class), any())).thenReturn(testUser);
    when(imputationService.queryForWorkspaceApps()).thenReturn(null);
  }

  @Test
  void testGetPipelines() throws Exception {
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
  void getPipeline() throws Exception {
    String pipelineIdString = TestUtils.TEST_PIPELINE_1.getPipelineId();
    PipelinesEnum pipelineIdEnum = PipelinesEnum.IMPUTATION;

    when(pipelinesServiceMock.getPipeline(pipelineIdEnum)).thenReturn(TestUtils.TEST_PIPELINE_1);

    MvcResult result =
        mockMvc
            .perform(get("/api/pipelines/v1alpha1/" + pipelineIdString))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiPipeline response =
        new ObjectMapper().readValue(result.getResponse().getContentAsString(), ApiPipeline.class);

    assertEquals(pipelineIdString, response.getPipelineId());
  }

  @Test
  void getPipelineCaseInsensitive() throws Exception {
    String pipelineIdString = "iMpUtAtIoN";
    PipelinesEnum pipelineIdEnum = PipelinesEnum.IMPUTATION;

    when(pipelinesServiceMock.getPipeline(pipelineIdEnum)).thenReturn(TestUtils.TEST_PIPELINE_1);

    MvcResult result =
        mockMvc
            .perform(get("/api/pipelines/v1alpha1/" + pipelineIdString))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiPipeline response =
        new ObjectMapper().readValue(result.getResponse().getContentAsString(), ApiPipeline.class);

    assertEquals(pipelineIdEnum.getValue(), response.getPipelineId());
  }

  @Test
  void getPipeline_badPipeline() throws Exception {
    String pipelineIdString = "bad-pipeline-id";

    mockMvc
        .perform(get("/api/pipelines/v1alpha1/" + pipelineIdString))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(InvalidPipelineException.class, result.getResolvedException()));
  }

  @Test
  void testCreateJobImputationPipelineRunning() throws Exception {
    String pipelineIdString = PipelinesEnum.IMPUTATION.getValue();
    String description = "description for testCreateJobImputationPipelineRunning";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);

    // the mocks
    when(imputationService.createImputationJob(
            testUser.getSubjectId(), description, testPipelineVersion, testPipelineInputs))
        .thenReturn(jobId);
    when(jobServiceMock.retrieveJob(jobId, testUser.getSubjectId()))
        .thenReturn(
            StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.RUNNING, jobId));

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelines/v1alpha1/%s", pipelineIdString))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isAccepted())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(ApiJobReport.StatusEnum.RUNNING, response.getJobReport().getStatus());
  }

  @Test
  void testCreateJobImputationPipelineNullDescriptionOk() throws Exception {
    String pipelineIdString = PipelinesEnum.IMPUTATION.getValue();
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), null);

    // the mocks
    when(imputationService.createImputationJob(
            testUser.getSubjectId(), null, testPipelineVersion, testPipelineInputs))
        .thenReturn(jobId);
    when(jobServiceMock.retrieveJob(jobId, testUser.getSubjectId()))
        .thenReturn(
            StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.RUNNING, jobId));

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelines/v1alpha1/%s", pipelineIdString))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isAccepted())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(ApiJobReport.StatusEnum.RUNNING, response.getJobReport().getStatus());
  }

  @Test
  void testCreateJobImputationPipelineCompletedSuccess() throws Exception {
    String pipelineIdString = PipelinesEnum.IMPUTATION.getValue();
    String description = "description for testCreateJobImputationPipelineCompletedSuccess";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);

    // the mocks
    when(imputationService.createImputationJob(
            testUser.getSubjectId(), description, testPipelineVersion, testPipelineInputs))
        .thenReturn(jobId);
    when(jobServiceMock.retrieveJob(jobId, testUser.getSubjectId()))
        .thenReturn(
            StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.SUCCESS, jobId));

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelines/v1alpha1/%s", pipelineIdString))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(ApiJobReport.StatusEnum.SUCCEEDED, response.getJobReport().getStatus());
  }

  private String createTestJobPostBody(String jobId, String description)
      throws JsonProcessingException {
    ApiJobControl apiJobControl = new ApiJobControl();
    apiJobControl.setId(jobId);
    ApiCreateJobRequestBody postBody =
        new ApiCreateJobRequestBody()
            //             .jobControl(apiJobControl)
            .pipelineVersion(testPipelineVersion)
            .pipelineInputs(testPipelineInputs)
            .description(description);
    return MockMvcUtils.convertToJsonString(postBody);
  }

  @Test
  void testCreateJobImputationPipelineCaseInsensitive() throws Exception {
    String pipelineIdString = "iMpUtAtIoN";
    String description = "description for testCreateJobImputationPipelineCaseInsensitive";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);

    // the mocks
    when(imputationService.createImputationJob(
            testUser.getSubjectId(), description, testPipelineVersion, testPipelineInputs))
        .thenReturn(jobId);
    when(jobServiceMock.retrieveJob(jobId, testUser.getSubjectId()))
        .thenReturn(
            StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.SUCCESS, jobId));

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelines/v1alpha1/%s", pipelineIdString))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
  }

  @Test
  void testCreateJobBadPipeline() throws Exception {
    String pipelineIdString = "bad-pipeline-id";
    String description = "description for testCreateJobBadPipeline";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);

    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1alpha1/%s", pipelineIdString))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(InvalidPipelineException.class, result.getResolvedException()));
  }

  @Test
  void testCreateJobMissingJobId() throws Exception {
    String pipelineIdString = PipelinesEnum.IMPUTATION.getValue();
    ApiCreateJobRequestBody postBody =
        new ApiCreateJobRequestBody()
            .pipelineVersion(testPipelineVersion)
            .pipelineInputs(testPipelineInputs)
            .description("description for testCreateJobMissingJobId");
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1alpha1/%s", pipelineIdString))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(
                    MissingRequiredFieldException.class, result.getResolvedException()));
  }

  @Test
  void testCreateImputationJobStairwayError() throws Exception {
    String pipelineIdString = PipelinesEnum.IMPUTATION.getValue();
    String description = "description for testCreateImputationJobStairwayError";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);

    // the mocks - one error that can happen is a MissingRequiredFieldException from Stairway
    when(imputationService.createImputationJob(
            testUser.getSubjectId(), description, testPipelineVersion, testPipelineInputs))
        .thenThrow(new InternalStairwayException("some message"));

    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1alpha1/%s", pipelineIdString))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isInternalServerError())
        .andExpect(
            result ->
                assertInstanceOf(InternalStairwayException.class, result.getResolvedException()));
  }

  @Test
  void testGetPipelineJobs() throws Exception {
    String pipelineIdString = "imputation";
    PipelinesEnum pipelineIdEnum = PipelinesEnum.IMPUTATION;

    UUID jobId1 = UUID.randomUUID();
    UUID jobId2 = UUID.randomUUID();
    UUID jobId3 = UUID.randomUUID();
    EnumeratedJob job1Running =
        new EnumeratedJob()
            .flightState(
                StairwayTestUtils.constructFlightStateWithStatusAndId(
                    FlightStatus.RUNNING, jobId1));
    EnumeratedJob job2Success =
        new EnumeratedJob()
            .flightState(
                StairwayTestUtils.constructFlightStateWithStatusAndId(
                    FlightStatus.SUCCESS, jobId2));
    FlightState flightStateError =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.ERROR, jobId3);
    flightStateError.setException(new Exception("Test exception"));
    EnumeratedJob job3Error = new EnumeratedJob().flightState(flightStateError);

    EnumeratedJobs allJobs =
        new EnumeratedJobs().results(List.of(job1Running, job2Success, job3Error)).totalResults(3);

    when(jobServiceMock.enumerateJobs(testUser.getSubjectId(), 10, null, pipelineIdEnum))
        .thenReturn(allJobs);

    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/pipelines/v1alpha1/%s/jobs", pipelineIdString)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiGetJobsResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiGetJobsResponse.class);

    assertEquals(3, response.getTotalResults());
    assertEquals(3, response.getResults().size());
    assertArrayEquals(
        new String[] {jobId1.toString(), jobId2.toString(), jobId3.toString()},
        response.getResults().stream().map(ApiJobReport::getId).toArray());
    assertArrayEquals(
        new ApiJobReport.StatusEnum[] {
          ApiJobReport.StatusEnum.RUNNING,
          ApiJobReport.StatusEnum.SUCCEEDED,
          ApiJobReport.StatusEnum.FAILED
        },
        response.getResults().stream().map(ApiJobReport::getStatus).toArray());
  }

  @Test
  void testGetPipelineJobs_badPipeline() throws Exception {
    String pipelineIdString = "bad-pipeline-id";

    mockMvc
        .perform(get(String.format("/api/pipelines/v1alpha1/%s/jobs", pipelineIdString)))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(InvalidPipelineException.class, result.getResolvedException()));
  }
}
