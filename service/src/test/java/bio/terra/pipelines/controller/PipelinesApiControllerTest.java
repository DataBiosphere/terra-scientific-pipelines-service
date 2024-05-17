package bio.terra.pipelines.controller;

import static bio.terra.pipelines.testutils.MockMvcUtils.getTestPipeline;
import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.*;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.post;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.common.exception.ValidationException;
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
import bio.terra.pipelines.db.entities.PipelineInputDefinition;
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
import jakarta.servlet.http.HttpServletRequest;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.http.MediaType;
import org.springframework.http.converter.HttpMessageNotReadableException;
import org.springframework.mock.web.MockHttpServletRequest;
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
  @MockBean IngressConfiguration ingressConfiguration;

  @Autowired private MockMvc mockMvc;

  private final List<Pipeline> testPipelineList =
      List.of(TestUtils.TEST_PIPELINE_1, TestUtils.TEST_PIPELINE_2);
  private final SamUser testUser = MockMvcUtils.TEST_SAM_USER;
  private final String testPipelineVersion = TestUtils.TEST_PIPELINE_VERSION_1;
  private final Pipeline testPipeline = TestUtils.TEST_PIPELINE_1;
  private final Map<String, Object> testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID newJobId = TestUtils.TEST_NEW_UUID;
  private final String fullResultURL = TestUtils.TEST_RESULT_URL;

  @BeforeEach
  void beforeEach() {
    when(ingressConfiguration.getDomainName()).thenReturn(TestUtils.TEST_DOMAIN);
    when(samUserFactoryMock.from(any(HttpServletRequest.class), any())).thenReturn(testUser);
    when(pipelinesServiceMock.getPipeline(any())).thenReturn(getTestPipeline());
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

    assertEquals(testPipelineList.size(), response.size());
  }

  @Test
  void getPipelineOk() throws Exception {
    String pipelineName = TestUtils.TEST_PIPELINE_1.getName();
    PipelinesEnum pipelineNameEnum = PipelinesEnum.IMPUTATION_BEAGLE;

    when(pipelinesServiceMock.getPipeline(pipelineNameEnum)).thenReturn(TestUtils.TEST_PIPELINE_1);

    MvcResult result =
        mockMvc
            .perform(get("/api/pipelines/v1/" + pipelineName))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiPipelineWithDetails response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiPipelineWithDetails.class);

    assertEquals(pipelineName, response.getPipelineName());

    // check that the response only includes user-provided inputs
    assertEquals(
        TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST.stream()
            .filter(PipelineInputDefinition::getUserProvided)
            .toList()
            .size(),
        response.getInputs().size());
    for (ApiPipelineUserProvidedInputDefinition p : response.getInputs()) {
      // find the matching input definition in test pipeline inputs list and check if it's user
      // provided
      assertTrue(
          TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST.stream()
              .anyMatch(i -> i.getName().equals(p.getName()) && i.getUserProvided().equals(true)));
    }
  }

  @Test
  void getPipelineCaseInsensitive() throws Exception {
    String pipelineName = "ImpuTatioN_bEaGlE";
    PipelinesEnum pipelineNameEnum = PipelinesEnum.IMPUTATION_BEAGLE;

    when(pipelinesServiceMock.getPipeline(pipelineNameEnum)).thenReturn(TestUtils.TEST_PIPELINE_1);

    MvcResult result =
        mockMvc
            .perform(get("/api/pipelines/v1/" + pipelineName))
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
        .perform(get("/api/pipelines/v1/" + pipelineName))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(InvalidPipelineException.class, result.getResolvedException()));
  }

  // createPipelineJob tests

  @Test
  void createJobImputationPipelineRunning() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String description = "description for testCreateJobImputationPipelineRunning";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);
    FlightMap inputParameters = StairwayTestUtils.constructCreateJobInputs(new FlightMap());
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.RUNNING, jobId, inputParameters, new FlightMap());
    String resultPath = "https://" + TestUtils.TEST_DOMAIN + "/result/" + jobId;

    // the mocks
    when(imputationService.createImputationJob(
            jobId,
            testUser.getSubjectId(),
            description,
            testPipeline,
            testPipelineInputs,
            resultPath))
        .thenReturn(jobId);
    when(jobServiceMock.retrieveJob(
            jobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_BEAGLE))
        .thenReturn(flightState);

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelines/v1/%s", pipelineName))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isAccepted())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(fullResultURL, response.getJobReport().getResultURL());
    assertEquals(ApiJobReport.StatusEnum.RUNNING, response.getJobReport().getStatus());
  }

  @Test
  void createJobImputationPipelineNoDescriptionOk() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), "");
    String resultPath = "https://" + TestUtils.TEST_DOMAIN + "/result/" + jobId;

    // the mocks
    when(imputationService.createImputationJob(
            jobId, testUser.getSubjectId(), "", testPipeline, testPipelineInputs, resultPath))
        .thenReturn(jobId);
    when(jobServiceMock.retrieveJob(
            jobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_BEAGLE))
        .thenReturn(
            StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.RUNNING, jobId));

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelines/v1/%s", pipelineName))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isAccepted())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(fullResultURL, response.getJobReport().getResultURL());
    assertEquals(ApiJobReport.StatusEnum.RUNNING, response.getJobReport().getStatus());
  }

  @Test
  void createJobImputationPipelineCompletedSuccess() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String description = "description for testCreateJobImputationPipelineCompletedSuccess";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);
    String resultPath = "https://" + TestUtils.TEST_DOMAIN + "/result/" + jobId;

    // the mocks
    when(imputationService.createImputationJob(
            jobId,
            testUser.getSubjectId(),
            description,
            testPipeline,
            testPipelineInputs,
            resultPath))
        .thenReturn(jobId);
    when(jobServiceMock.retrieveJob(
            jobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_BEAGLE))
        .thenReturn(
            StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.SUCCESS, jobId));

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelines/v1/%s", pipelineName))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(fullResultURL, response.getJobReport().getResultURL());
    assertEquals(ApiJobReport.StatusEnum.SUCCEEDED, response.getJobReport().getStatus());
  }

  @Test
  void createJobImputationPipelineCaseInsensitive() throws Exception {
    String pipelineName = "iMpUtAtIoN_BEAGle";
    String description = "description for testCreateJobImputationPipelineCaseInsensitive";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);
    String resultPath = "https://" + TestUtils.TEST_DOMAIN + "/result/" + jobId;

    // the mocks
    when(imputationService.createImputationJob(
            jobId,
            testUser.getSubjectId(),
            description,
            testPipeline,
            testPipelineInputs,
            resultPath))
        .thenReturn(jobId);
    when(jobServiceMock.retrieveJob(
            jobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_BEAGLE))
        .thenReturn(
            StairwayTestUtils.constructFlightStateWithStatusAndId(FlightStatus.SUCCESS, jobId));

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelines/v1/%s", pipelineName))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(fullResultURL, response.getJobReport().getResultURL());
  }

  @Test
  void createJobBadPipeline() throws Exception {
    String pipelineName = "bad-pipeline-id";
    String description = "description for testCreateJobBadPipeline";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);

    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1/%s", pipelineName))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(InvalidPipelineException.class, result.getResolvedException()));
  }

  @Test
  void createJobMissingJobControl() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    ApiPipelineUserProvidedInputs userProvidedInputs = new ApiPipelineUserProvidedInputs();
    userProvidedInputs.putAll(testPipelineInputs);
    ApiCreateJobRequestBody postBody =
        new ApiCreateJobRequestBody()
            .pipelineVersion(testPipelineVersion)
            .pipelineInputs(userProvidedInputs)
            .description("description for testCreateJobMissingJobId");
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // Spring will catch the missing jobControl and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1/%s", pipelineName))
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
  void createJobMissingJobId() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    ApiJobControl apiJobControl = new ApiJobControl();
    ApiPipelineUserProvidedInputs userProvidedInputs = new ApiPipelineUserProvidedInputs();
    userProvidedInputs.putAll(testPipelineInputs);
    ApiCreateJobRequestBody postBody =
        new ApiCreateJobRequestBody()
            .jobControl(apiJobControl)
            .pipelineVersion(testPipelineVersion)
            .pipelineInputs(userProvidedInputs)
            .description("description for testCreateJobMissingJobId");
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // Spring will catch the missing job id and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1/%s", pipelineName))
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
  void createJobMissingMultipleRequiredFields() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String stringifiedInputs = MockMvcUtils.convertToJsonString(testPipelineInputs);
    String postBodyAsJson =
        String.format(
            "{\"pipelineInputs\":%s,\"description\":\"test description for testCreateJobMissingMultipleRequiredFields\"}",
            stringifiedInputs);

    // Spring will catch the missing fields and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1/%s", pipelineName))
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
  void createJobBadPipelineInputs() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String description = "description for testCreateJobBadPipelineInputs";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);

    // the mocks
    doThrow(new ValidationException("some message"))
        .when(pipelinesServiceMock)
        .validateUserProvidedInputs(PipelinesEnum.IMPUTATION_BEAGLE, testPipelineInputs);

    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1/%s", pipelineName))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isBadRequest())
        .andExpect(
            result -> assertInstanceOf(ValidationException.class, result.getResolvedException()));
  }

  @Test
  void createJobBadJobId() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String postBodyAsJson =
        createTestJobPostBody("this-is-not-a-uuid", "description for testCreateJobMissingJobId");

    // Spring will catch the non-uuid jobId and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1/%s", pipelineName))
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
  void createImputationJobStairwayError() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String description = "description for testCreateImputationJobStairwayError";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestJobPostBody(jobId.toString(), description);
    String resultPath = "https://" + TestUtils.TEST_DOMAIN + "/result/" + jobId;

    // the mocks - one error that can happen is a MissingRequiredFieldException from Stairway
    when(imputationService.createImputationJob(
            jobId,
            testUser.getSubjectId(),
            description,
            getTestPipeline(),
            testPipelineInputs,
            resultPath))
        .thenThrow(new InternalStairwayException("some message"));

    mockMvc
        .perform(
            post(String.format("/api/pipelines/v1/%s", pipelineName))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isInternalServerError())
        .andExpect(
            result ->
                assertInstanceOf(InternalStairwayException.class, result.getResolvedException()));
  }

  // getPipelineJobs tests

  @Test
  void getPipelineJobsOk() throws Exception {
    String pipelineName = "imputation_beagle";
    PipelinesEnum pipelineNameEnum = PipelinesEnum.IMPUTATION_BEAGLE;

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
            .perform(get(String.format("/api/pipelines/v1/%s/jobs", pipelineName)))
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
  void getPipelineJobsBadPipeline() throws Exception {
    String badPipelineName = "bad-pipeline-id";

    mockMvc
        .perform(get(String.format("/api/pipelines/v1/%s/jobs", badPipelineName)))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(InvalidPipelineException.class, result.getResolvedException()));
  }

  // getPipelineJobResult tests

  @Test
  void getPipelineJobResultDoneSuccess() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String jobIdString = newJobId.toString();
    String jobResultValue = "job result value";

    JobApiUtils.AsyncJobResult<String> jobResult =
        new JobApiUtils.AsyncJobResult<String>()
            .jobReport(
                new ApiJobReport()
                    .id(newJobId.toString())
                    .status(ApiJobReport.StatusEnum.SUCCEEDED))
            .result(jobResultValue);

    // the mocks
    when(jobServiceMock.retrieveAsyncJobResult(
            newJobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_BEAGLE, String.class, null))
        .thenReturn(jobResult);

    MvcResult result =
        mockMvc
            .perform(
                get(String.format("/api/pipelines/v1/%s/result/%s", pipelineName, jobIdString)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);

    // response should include the job report and pipeline output object
    assertEquals(newJobId.toString(), response.getJobReport().getId());
    assertEquals(jobResultValue, response.getPipelineOutput());
    assertNull(response.getErrorReport());
  }

  @Test
  void getPipelineJobResultDoneFailed() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String jobIdString = newJobId.toString();
    String errorMessage = "test exception message";
    Integer statusCode = 500;

    ApiErrorReport errorReport = new ApiErrorReport().message(errorMessage).statusCode(statusCode);

    JobApiUtils.AsyncJobResult<String> jobResult =
        new JobApiUtils.AsyncJobResult<String>()
            .jobReport(
                new ApiJobReport()
                    .id(newJobId.toString())
                    .status(ApiJobReport.StatusEnum.FAILED)
                    .statusCode(statusCode))
            .errorReport(errorReport);

    // the mocks
    when(jobServiceMock.retrieveAsyncJobResult(
            newJobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_BEAGLE, String.class, null))
        .thenReturn(jobResult);

    MvcResult result =
        mockMvc
            .perform(
                get(String.format("/api/pipelines/v1/%s/result/%s", pipelineName, jobIdString)))
            .andExpect(status().isOk()) // the call itself should return a 200
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);

    // response should include the error report and no pipeline output object
    assertEquals(newJobId.toString(), response.getJobReport().getId());
    assertNull(response.getPipelineOutput());
    assertEquals(statusCode, response.getJobReport().getStatusCode());
    assertEquals(errorMessage, response.getErrorReport().getMessage());
  }

  @Test
  void getPipelineJobResultRunning() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String jobIdString = newJobId.toString();
    Integer statusCode = 202;
    JobApiUtils.AsyncJobResult<String> jobResult =
        new JobApiUtils.AsyncJobResult<String>()
            .jobReport(
                new ApiJobReport()
                    .id(newJobId.toString())
                    .status(ApiJobReport.StatusEnum.RUNNING)
                    .statusCode(statusCode));

    // the mocks
    when(jobServiceMock.retrieveAsyncJobResult(
            newJobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_BEAGLE, String.class, null))
        .thenReturn(jobResult);

    MvcResult result =
        mockMvc
            .perform(
                get(String.format("/api/pipelines/v1/%s/result/%s", pipelineName, jobIdString)))
            .andExpect(status().isAccepted())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiCreateJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiCreateJobResponse.class);

    // response should include the job report and no error report or pipeline output object
    assertEquals(newJobId.toString(), response.getJobReport().getId());
    assertEquals(statusCode, response.getJobReport().getStatusCode());
    assertNull(response.getPipelineOutput());
    assertNull(response.getErrorReport());
  }

  @Test
  void getAsyncResultEndpointHttps() {
    String testServletPath = "test/path";

    MockHttpServletRequest request = new MockHttpServletRequest();
    request.setServletPath(testServletPath);

    UUID jobId = newJobId;
    // the function prepends https:// and the domain to the path, and append "result" and the jobId
    String expectedResultEndpoint =
        String.format("https://%s/%s/result/%s", TestUtils.TEST_DOMAIN, testServletPath, jobId);

    assertEquals(
        expectedResultEndpoint,
        PipelinesApiController.getAsyncResultEndpoint(ingressConfiguration, request, jobId));
  }

  @Test
  void getAsyncResultEndpointHttp() {
    String testServletPath = "test/path";

    MockHttpServletRequest request = new MockHttpServletRequest();
    request.setServletPath(testServletPath);

    // override this mock to return localhost
    String localhostDomain = "localhost:8080";
    when(ingressConfiguration.getDomainName()).thenReturn(localhostDomain);

    UUID jobId = newJobId;
    // for localhost, the function prepends http:// and the domain to the path, and append "result"
    // and the jobId
    String expectedResultEndpoint =
        String.format("http://%s/%s/result/%s", localhostDomain, testServletPath, jobId);

    assertEquals(
        expectedResultEndpoint,
        PipelinesApiController.getAsyncResultEndpoint(ingressConfiguration, request, jobId));
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
