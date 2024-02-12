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
import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.app.controller.GlobalExceptionHandler;
import bio.terra.pipelines.app.controller.JobApiUtils;
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
import bio.terra.stairway.FlightMap;
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
import org.springframework.boot.test.mock.mockito.SpyBean;
import org.springframework.http.MediaType;
import org.springframework.http.converter.HttpMessageNotReadableException;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.MvcResult;
import org.springframework.test.web.servlet.result.MockMvcResultMatchers;
import org.springframework.web.bind.MethodArgumentNotValidException;

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
  @SpyBean JobApiUtils jobApiUtils;
  @MockBean IngressConfiguration ingressConfiguration;

  @Autowired private MockMvc mockMvc;

  private final List<Pipeline> testPipelineList =
      List.of(TestUtils.TEST_PIPELINE_1, TestUtils.TEST_PIPELINE_2);
  private final SamUser testUser = MockMvcUtils.TEST_SAM_USER;
  private final String testPipelineVersion = TestUtils.TEST_PIPELINE_VERSION_1;
  private final Pipeline testPipeline = TestUtils.TEST_PIPELINE_1;
  private final Object testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID newJobId = TestUtils.TEST_NEW_UUID;
  private final String testResultPath = TestUtils.TEST_RESULT_PATH;
  private final String testResultPathURL = String.format("http://localhost/%s", testResultPath);

  @BeforeEach
  void beforeEach() {
    jobApiUtils = spy(new JobApiUtils(jobServiceMock, ingressConfiguration));
    when(ingressConfiguration.getDomainName()).thenReturn("localhost");
    when(samUserFactoryMock.from(any(HttpServletRequest.class), any())).thenReturn(testUser);
    when(imputationService.queryForWorkspaceApps(any())).thenReturn(null);
    when(pipelinesServiceMock.getPipeline(any())).thenReturn(getTestPipeline());
  }

  // getPipeline tests

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
    String pipelineName = TestUtils.TEST_PIPELINE_1.getName();
    PipelinesEnum pipelineNameEnum = PipelinesEnum.IMPUTATION_MINIMAC4;

    when(pipelinesServiceMock.getPipeline(pipelineNameEnum)).thenReturn(TestUtils.TEST_PIPELINE_1);

    MvcResult result =
        mockMvc
            .perform(get("/api/pipelines/v1alpha1/" + pipelineName))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiPipeline response =
        new ObjectMapper().readValue(result.getResponse().getContentAsString(), ApiPipeline.class);

    assertEquals(pipelineName, response.getPipelineName());
  }

  @Test
  void getPipelineCaseInsensitive() throws Exception {
    String pipelineName = "ImpuTatioN_MIniMac4";
    PipelinesEnum pipelineNameEnum = PipelinesEnum.IMPUTATION_MINIMAC4;

    when(pipelinesServiceMock.getPipeline(pipelineNameEnum)).thenReturn(TestUtils.TEST_PIPELINE_1);

    MvcResult result =
        mockMvc
            .perform(get("/api/pipelines/v1alpha1/" + pipelineName))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiPipeline response =
        new ObjectMapper().readValue(result.getResponse().getContentAsString(), ApiPipeline.class);

    assertEquals(pipelineNameEnum.getValue(), response.getPipelineName());
  }

  @Test
  void getPipeline_badPipeline() throws Exception {
    String pipelineName = "bad-pipeline-id";

    mockMvc
        .perform(get("/api/pipelines/v1alpha1/" + pipelineName))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(InvalidPipelineException.class, result.getResolvedException()));
  }

  // createPipelineJob tests

  @Test
  void testCreateJobImputationPipelineRunning() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_MINIMAC4.getValue();
    String description = "description for testCreateJobImputationPipelineRunning";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);
    FlightMap inputParameters = StairwayTestUtils.constructCreateJobInputs(new FlightMap());
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.RUNNING, jobId, inputParameters, new FlightMap());

    // the mocks
    when(imputationService.createImputationJob(
            jobId,
            testUser.getSubjectId(),
            description,
            testPipeline,
            testPipelineInputs,
            testResultPath))
        .thenReturn(jobId);
    when(jobServiceMock.retrieveJob(
            jobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_MINIMAC4))
        .thenReturn(flightState);

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelines/v1alpha1/%s", pipelineName))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isAccepted())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(testResultPathURL, response.getJobReport().getResultURL());
    assertEquals(ApiJobReport.StatusEnum.RUNNING, response.getJobReport().getStatus());
  }

  @Test
  void testCreateJobImputationPipelineNoDescriptionOk() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_MINIMAC4.getValue();
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), "");

    // the mocks
    when(imputationService.createImputationJob(
            jobId, testUser.getSubjectId(), "", testPipeline, testPipelineInputs, testResultPath))
        .thenReturn(jobId);
    when(jobServiceMock.retrieveJob(
            jobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_MINIMAC4))
        .thenReturn(
            StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.RUNNING, jobId));

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelines/v1alpha1/%s", pipelineName))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isAccepted())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(testResultPathURL, response.getJobReport().getResultURL());
    assertEquals(ApiJobReport.StatusEnum.RUNNING, response.getJobReport().getStatus());
  }

  @Test
  void testCreateJobImputationPipelineCompletedSuccess() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_MINIMAC4.getValue();
    String description = "description for testCreateJobImputationPipelineCompletedSuccess";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);

    // the mocks
    when(imputationService.createImputationJob(
            jobId,
            testUser.getSubjectId(),
            description,
            testPipeline,
            testPipelineInputs,
            testResultPath))
        .thenReturn(jobId);
    when(jobServiceMock.retrieveJob(
            jobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_MINIMAC4))
        .thenReturn(
            StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.SUCCESS, jobId));

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelines/v1alpha1/%s", pipelineName))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(testResultPathURL, response.getJobReport().getResultURL());
    assertEquals(ApiJobReport.StatusEnum.SUCCEEDED, response.getJobReport().getStatus());
  }

  @Test
  void testCreateJobImputationPipelineCaseInsensitive() throws Exception {
    String pipelineName = "iMpUtAtIoN_MINImac4";
    String description = "description for testCreateJobImputationPipelineCaseInsensitive";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);

    // the mocks
    when(imputationService.createImputationJob(
            jobId,
            testUser.getSubjectId(),
            description,
            testPipeline,
            testPipelineInputs,
            testResultPath))
        .thenReturn(jobId);
    when(jobServiceMock.retrieveJob(
            jobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_MINIMAC4))
        .thenReturn(
            StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.SUCCESS, jobId));

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelines/v1alpha1/%s", pipelineName))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(testResultPathURL, response.getJobReport().getResultURL());
  }

  @Test
  void testCreateJobBadPipeline() throws Exception {
    String pipelineName = "bad-pipeline-id";
    String description = "description for testCreateJobBadPipeline";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);

    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1alpha1/%s", pipelineName))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(InvalidPipelineException.class, result.getResolvedException()));
  }

  @Test
  void testCreateJobMissingJobControl() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_MINIMAC4.getValue();
    ApiCreateJobRequestBody postBody =
        new ApiCreateJobRequestBody()
            .pipelineVersion(testPipelineVersion)
            .pipelineInputs(testPipelineInputs)
            .description("description for testCreateJobMissingJobId");
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // Spring will catch the missing jobControl and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1alpha1/%s", pipelineName))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(
                    MethodArgumentNotValidException.class, result.getResolvedException()))
        .andExpect(
            MockMvcResultMatchers.jsonPath("$.message")
                .value("Request could not be parsed or was invalid: jobControl must not be null"));
  }

  @Test
  void testCreateJobMissingJobId() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_MINIMAC4.getValue();
    ApiJobControl apiJobControl = new ApiJobControl();
    ApiCreateJobRequestBody postBody =
        new ApiCreateJobRequestBody()
            .jobControl(apiJobControl)
            .pipelineVersion(testPipelineVersion)
            .pipelineInputs(testPipelineInputs)
            .description("description for testCreateJobMissingJobId");
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // Spring will catch the missing job id and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1alpha1/%s", pipelineName))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(
                    MethodArgumentNotValidException.class, result.getResolvedException()))
        .andExpect(
            MockMvcResultMatchers.jsonPath("$.message")
                .value(
                    "Request could not be parsed or was invalid: jobControl.id must not be null"));
  }

  @Test
  void testCreateJobMissingMultipleRequiredFields() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_MINIMAC4.getValue();
    String stringifiedInputs = MockMvcUtils.convertToJsonString(testPipelineInputs);
    String postBodyAsJson =
        String.format(
            "{\"pipelineInputs\":%s,\"description\":\"test description for testCreateJobMissingMultipleRequiredFields\"}",
            stringifiedInputs);

    // Spring will catch the missing fields and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1alpha1/%s", pipelineName))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(
                    MethodArgumentNotValidException.class, result.getResolvedException()))
        .andExpect(
            MockMvcResultMatchers.jsonPath("$.message")
                .value(
                    "Request could not be parsed or was invalid: jobControl must not be null; "
                        + "pipelineVersion must not be null"));
  }

  @Test
  void testCreateJobBadJobId() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_MINIMAC4.getValue();
    String postBodyAsJson =
        createTestJobPostBody("this-is-not-a-uuid", "description for testCreateJobMissingJobId");

    // Spring will catch the non-uuid jobId and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1alpha1/%s", pipelineName))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(
                    HttpMessageNotReadableException.class, result.getResolvedException()))
        .andExpect(
            MockMvcResultMatchers.jsonPath("$.message")
                .value(
                    "JSON parse error: Cannot deserialize value of type `java.util.UUID` from String \"this-is-not-a-uuid\": UUID has to be represented by standard 36-char representation"));
  }

  @Test
  void testCreateImputationJobStairwayError() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_MINIMAC4.getValue();
    String description = "description for testCreateImputationJobStairwayError";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);
    String resultPath = "/result/" + jobId;

    // the mocks - one error that can happen is a MissingRequiredFieldException from Stairway
    when(imputationService.createImputationJob(
            jobId,
            testUser.getSubjectId(),
            description,
            MockMvcUtils.getTestPipeline(),
            testPipelineInputs,
            resultPath))
        .thenThrow(new InternalStairwayException("some message"));

    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1alpha1/%s", pipelineName))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isInternalServerError())
        .andExpect(
            result ->
                assertInstanceOf(InternalStairwayException.class, result.getResolvedException()));
  }

  // getPipelineJobs tests

  @Test
  void testGetPipelineJobs() throws Exception {
    String pipelineName = "imputation_minimac4";
    PipelinesEnum pipelineNameEnum = PipelinesEnum.IMPUTATION_MINIMAC4;

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

    when(jobServiceMock.enumerateJobs(testUser.getSubjectId(), 10, null, pipelineNameEnum))
        .thenReturn(allJobs);

    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/pipelines/v1alpha1/%s/jobs", pipelineName)))
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
    String badPipelineName = "bad-pipeline-id";

    mockMvc
        .perform(get(String.format("/api/pipelines/v1alpha1/%s/jobs", badPipelineName)))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(InvalidPipelineException.class, result.getResolvedException()));
  }

  // getPipelineJobResult tests

  @Test
  void testGetPipelineJobResultDoneSuccess() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_MINIMAC4.getValue();
    String jobIdString = newJobId.toString();
    String jobResultValue = "job result value";
    FlightState expectedFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.SUCCESS, newJobId);

    JobService.JobResultOrException<String> resultOrException =
        new JobService.JobResultOrException<String>().result(jobResultValue);

    // the mocks
    when(jobServiceMock.retrieveJob(
            newJobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_MINIMAC4))
        .thenReturn(expectedFlightState);
    when(jobServiceMock.retrieveJobResult(newJobId, String.class, null))
        .thenReturn(resultOrException);

    MvcResult result =
        mockMvc
            .perform(
                get(
                    String.format(
                        "/api/pipelines/v1alpha1/%s/result/%s", pipelineName, jobIdString)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);

    // response should include the job report and pipeline output object
    assertEquals(newJobId.toString(), response.getJobReport().getId());
    assertNotNull(response.getPipelineOutput());
    assertNull(response.getErrorReport());
  }

  @Test
  void testGetPipelineJobResultDoneFailed() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_MINIMAC4.getValue();
    String jobIdString = newJobId.toString();

    FlightState expectedFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.ERROR, newJobId);
    expectedFlightState.setException(new Exception("Test exception"));

    JobService.JobResultOrException<String> resultOrException =
        new JobService.JobResultOrException<String>()
            .exception(new RuntimeException("Test exception"));

    // the mocks
    when(jobServiceMock.retrieveJob(
            newJobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_MINIMAC4))
        .thenReturn(expectedFlightState);
    when(jobServiceMock.retrieveJobResult(newJobId, String.class, null))
        .thenReturn(resultOrException);

    MvcResult result =
        mockMvc
            .perform(
                get(
                    String.format(
                        "/api/pipelines/v1alpha1/%s/result/%s", pipelineName, jobIdString)))
            .andExpect(status().isOk()) // the call itself should return a 200
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);

    // response should include the error report and no pipeline output object
    assertEquals(newJobId.toString(), response.getJobReport().getId());
    assertNull(response.getPipelineOutput());
    assertEquals(500, response.getJobReport().getStatusCode());
  }

  @Test
  void testGetPipelineJobResultRunning() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_MINIMAC4.getValue();
    String jobIdString = newJobId.toString();
    FlightState expectedFlightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.RUNNING, newJobId);

    // the mocks
    when(jobServiceMock.retrieveJob(
            newJobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_MINIMAC4))
        .thenReturn(expectedFlightState);

    MvcResult result =
        mockMvc
            .perform(
                get(
                    String.format(
                        "/api/pipelines/v1alpha1/%s/result/%s", pipelineName, jobIdString)))
            .andExpect(status().isAccepted())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);

    // response should include the job report and no error report or pipeline output object
    assertEquals(newJobId.toString(), response.getJobReport().getId());
    assertNull(response.getPipelineOutput());
    assertNull(response.getErrorReport());
  }

  // support methods

  private String createTestJobPostBody(String jobId, String description)
      throws JsonProcessingException {
    String stringifiedInputs = MockMvcUtils.convertToJsonString(testPipelineInputs);
    return String.format(
        "{\"jobControl\":{\"id\":\"%s\"},\"pipelineVersion\":\"%s\",\"pipelineInputs\":%s,\"description\":\"%s\"}",
        jobId, testPipelineVersion, stringifiedInputs, description);
  }
}
