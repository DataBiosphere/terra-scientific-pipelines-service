package bio.terra.pipelines.controller;

import static bio.terra.pipelines.testutils.MockMvcUtils.getTestPipeline;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.doNothing;
import static org.mockito.Mockito.doThrow;
import static org.mockito.Mockito.when;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.post;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.common.iam.BearerTokenFactory;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.app.controller.GlobalExceptionHandler;
import bio.terra.pipelines.app.controller.JobApiUtils;
import bio.terra.pipelines.app.controller.PipelineRunsApiController;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.dependencies.stairway.exception.InternalStairwayException;
import bio.terra.pipelines.generated.model.ApiAsyncPipelineRunResponse;
import bio.terra.pipelines.generated.model.ApiCreatePipelineRunRequestBody;
import bio.terra.pipelines.generated.model.ApiErrorReport;
import bio.terra.pipelines.generated.model.ApiJobControl;
import bio.terra.pipelines.generated.model.ApiJobReport;
import bio.terra.pipelines.generated.model.ApiPipelineRunOutput;
import bio.terra.pipelines.generated.model.ApiPipelineUserProvidedInputs;
import bio.terra.pipelines.service.PipelineRunsService;
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
import java.time.LocalDateTime;
import java.util.Map;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.http.HttpStatus;
import org.springframework.http.MediaType;
import org.springframework.http.converter.HttpMessageNotReadableException;
import org.springframework.mock.web.MockHttpServletRequest;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.MvcResult;
import org.springframework.test.web.servlet.result.MockMvcResultMatchers;
import org.springframework.web.bind.MethodArgumentNotValidException;

@ContextConfiguration(classes = {PipelineRunsApiController.class, GlobalExceptionHandler.class})
@WebMvcTest
class PipelineRunsApiControllerTest {
  @MockBean PipelinesService pipelinesServiceMock;
  @MockBean PipelineRunsService pipelineRunsServiceMock;
  @MockBean JobService jobServiceMock;
  @MockBean SamUserFactory samUserFactoryMock;
  @MockBean BearerTokenFactory bearerTokenFactory;
  @MockBean SamConfiguration samConfiguration;
  @MockBean IngressConfiguration ingressConfiguration;

  @Autowired private MockMvc mockMvc;

  private final SamUser testUser = MockMvcUtils.TEST_SAM_USER;
  private final String testPipelineVersion = TestUtils.TEST_PIPELINE_VERSION_1;
  private final Map<String, Object> testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID newJobId = TestUtils.TEST_NEW_UUID;
  private final LocalDateTime createdTime = LocalDateTime.now();
  private final LocalDateTime updatedTime = LocalDateTime.now();
  private final String testResultPath = TestUtils.TEST_RESULT_URL;
  private final Map<String, String> testOutput = TestUtils.TEST_PIPELINE_OUTPUTS;
  private final ObjectMapper objectMapper = new ObjectMapper();
  private final String testOutputString = objectMapper.writeValueAsString(testOutput);

  PipelineRunsApiControllerTest() throws JsonProcessingException {}

  @BeforeEach
  void beforeEach() {
    when(ingressConfiguration.getDomainName()).thenReturn(TestUtils.TEST_DOMAIN);
    when(samUserFactoryMock.from(any(HttpServletRequest.class), any())).thenReturn(testUser);
    when(pipelinesServiceMock.getPipeline(any())).thenReturn(getTestPipeline());
  }

  // createPipelineRun tests
  // createPipelineRun performs the following actions:
  // 1. extract user-provided information from the request
  // 2. validate the user-provided information
  // 3. write run info to pipeline_runs db and create a new job in Stairway (in a transaction)
  // 4. query the pipeline_runs table for the job and configure a response object

  @Test
  void createJobImputationPipelineRunning() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String description = "description for testCreateJobImputationPipelineRunning";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestPipelineRunPostBody(jobId.toString(), description);
    FlightMap inputParameters = StairwayTestUtils.constructCreateJobInputs(new FlightMap());
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.RUNNING, jobId, inputParameters, new FlightMap());
    PipelineRun testPipelineRun = createPipelineRunRunning();
    testPipelineRun.setDescription(description);

    ApiJobReport jobReport =
        new ApiJobReport()
            .id(flightState.getFlightId())
            .description(description)
            .status(ApiJobReport.StatusEnum.RUNNING)
            .statusCode(HttpStatus.ACCEPTED.value())
            .submitted(testPipelineRun.getCreated().toString())
            .completed(null)
            .resultURL(testPipelineRun.getResultUrl());

    // the mocks
    doNothing().when(pipelinesServiceMock).validateUserProvidedInputs(any(), any());
    when(pipelineRunsServiceMock.createPipelineRun(any(), any(), any(), any(), any(), any()))
        .thenReturn(testPipelineRun);
    when(jobServiceMock.retrieveJob(
            jobId, testUser.getSubjectId(), PipelinesEnum.IMPUTATION_BEAGLE))
        .thenReturn(flightState);
    when(jobServiceMock.retrieveAsyncJobResult(jobId, testUser.getSubjectId(), String.class, null))
        .thenReturn(new JobApiUtils.AsyncJobResult<String>().jobReport(jobReport).result(null));

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelineruns/v1/%s", pipelineName))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isAccepted())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAsyncPipelineRunResponse response =
        new ObjectMapper()
            .readValue(
                result.getResponse().getContentAsString(), ApiAsyncPipelineRunResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(testResultPath, response.getJobReport().getResultURL());
    assertEquals(ApiJobReport.StatusEnum.RUNNING, response.getJobReport().getStatus());
    assertEquals(createdTime.toString(), response.getJobReport().getSubmitted());
    assertNull(response.getPipelineOutput());
    assertNull(response.getJobReport().getCompleted());
  }

  @Test
  void createPipelineRunCompletedSuccessFromDb() throws Exception {
    // when the run is marked as is_success in the pipeline_runs table, we pull all info from that
    // table and not from stairway
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String description = "description for testCreateJobImputationPipelineCompletedSuccess";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestPipelineRunPostBody(jobId.toString(), description);
    PipelineRun pipelineRun = createPipelineRunCompleted(CommonPipelineRunStatusEnum.SUCCEEDED);
    ApiPipelineRunOutput apiPipelineRunOutput = new ApiPipelineRunOutput();
    apiPipelineRunOutput.putAll(testOutput);

    // the mocks
    doNothing().when(pipelinesServiceMock).validateUserProvidedInputs(any(), any());
    when(pipelineRunsServiceMock.createPipelineRun(any(), any(), any(), any(), any(), any()))
        .thenReturn(pipelineRun);
    when(pipelineRunsServiceMock.formatPipelineRunOutputs(pipelineRun))
        .thenReturn(apiPipelineRunOutput);

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelineruns/v1/%s", pipelineName))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAsyncPipelineRunResponse response =
        new ObjectMapper()
            .readValue(
                result.getResponse().getContentAsString(), ApiAsyncPipelineRunResponse.class);
    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(testResultPath, response.getJobReport().getResultURL());
    assertEquals(ApiJobReport.StatusEnum.SUCCEEDED, response.getJobReport().getStatus());
    assertEquals(createdTime.toString(), response.getJobReport().getSubmitted());
    assertEquals(updatedTime.toString(), response.getJobReport().getCompleted());
    assertEquals(testOutput, response.getPipelineOutput());
  }

  @Test
  void createPipelineRunMissingJobControl() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    ApiPipelineUserProvidedInputs userProvidedInputs = new ApiPipelineUserProvidedInputs();
    userProvidedInputs.putAll(testPipelineInputs);
    ApiCreatePipelineRunRequestBody postBody =
        new ApiCreatePipelineRunRequestBody()
            .pipelineVersion(testPipelineVersion)
            .pipelineInputs(userProvidedInputs)
            .description("description for testCreateJobMissingJobId");
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // Spring will catch the missing jobControl and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelineruns/v1/%s", pipelineName))
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
  void createPipelineRunMissingJobId() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    ApiJobControl apiJobControl = new ApiJobControl();
    ApiPipelineUserProvidedInputs userProvidedInputs = new ApiPipelineUserProvidedInputs();
    userProvidedInputs.putAll(testPipelineInputs);
    ApiCreatePipelineRunRequestBody postBody =
        new ApiCreatePipelineRunRequestBody()
            .jobControl(apiJobControl)
            .pipelineVersion(testPipelineVersion)
            .pipelineInputs(userProvidedInputs)
            .description("description for testCreateJobMissingJobId");
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // Spring will catch the missing job id and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelineruns/v1/%s", pipelineName))
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
  void createPipelineRunMissingMultipleRequiredFields() throws Exception {
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
            post(String.format("/api/pipelineruns/v1/%s", pipelineName))
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
  void createPipelineRunBadPipelineInputs() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String description = "description for testCreateJobBadPipelineInputs";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestPipelineRunPostBody(jobId.toString(), description);

    // the mocks
    doThrow(new ValidationException("some message"))
        .when(pipelinesServiceMock)
        .validateUserProvidedInputs(
            TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST, testPipelineInputs);

    mockMvc
        .perform(
            post(String.format("/api/pipelineruns/v1/%s", pipelineName))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isBadRequest())
        .andExpect(
            result -> assertInstanceOf(ValidationException.class, result.getResolvedException()));
  }

  @Test
  void createPipelineRunBadJobId() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String postBodyAsJson =
        createTestPipelineRunPostBody(
            "this-is-not-a-uuid", "description for testCreateJobMissingJobId");

    // Spring will catch the non-uuid jobId and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelineruns/v1/%s", pipelineName))
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
  void createPipelineRunDbError() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String description = "description for createPipelineRunDbError";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestPipelineRunPostBody(jobId.toString(), description);

    // the mocks
    doNothing().when(pipelinesServiceMock).validateUserProvidedInputs(any(), any());
    when(pipelineRunsServiceMock.createPipelineRun(
            getTestPipeline(),
            jobId,
            testUser.getSubjectId(),
            description,
            testPipelineInputs,
            testResultPath))
        .thenThrow(new RuntimeException("some message"));

    mockMvc
        .perform(
            post(String.format("/api/pipelineruns/v1/%s", pipelineName))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isInternalServerError())
        .andExpect(
            result -> assertInstanceOf(RuntimeException.class, result.getResolvedException()));
  }

  @Test
  void createImputationRunStairwayError() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String description = "description for testCreateImputationJobStairwayError";
    UUID jobId = newJobId;
    String postBodyAsJson = createTestPipelineRunPostBody(jobId.toString(), description);
    String resultPath = "https://" + TestUtils.TEST_DOMAIN + "/result/" + jobId;

    // the mocks - one error that can happen is a MissingRequiredFieldException from Stairway
    doNothing().when(pipelinesServiceMock).validateUserProvidedInputs(any(), any());
    when(pipelineRunsServiceMock.createPipelineRun(
            getTestPipeline(),
            jobId,
            testUser.getSubjectId(),
            description,
            testPipelineInputs,
            resultPath))
        .thenThrow(new InternalStairwayException("some message"));

    mockMvc
        .perform(
            post(String.format("/api/pipelineruns/v1/%s", pipelineName))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isInternalServerError())
        .andExpect(
            result ->
                assertInstanceOf(InternalStairwayException.class, result.getResolvedException()));
  }

  // getPipelineRunResult tests

  @Test
  void getPipelineRunResultDoneSuccess() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String jobIdString = newJobId.toString();
    PipelineRun pipelineRun = createPipelineRunCompleted(CommonPipelineRunStatusEnum.SUCCEEDED);
    ApiPipelineRunOutput apiPipelineRunOutput = new ApiPipelineRunOutput();
    apiPipelineRunOutput.putAll(testOutput);

    // the mocks
    doNothing().when(pipelinesServiceMock).validateUserProvidedInputs(any(), any());
    when(pipelineRunsServiceMock.getPipelineRun(newJobId, testUser.getSubjectId()))
        .thenReturn(pipelineRun);
    when(pipelineRunsServiceMock.formatPipelineRunOutputs(pipelineRun))
        .thenReturn(apiPipelineRunOutput);

    MvcResult result =
        mockMvc
            .perform(
                get(String.format("/api/pipelineruns/v1/%s/result/%s", pipelineName, jobIdString)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAsyncPipelineRunResponse response =
        new ObjectMapper()
            .readValue(
                result.getResponse().getContentAsString(), ApiAsyncPipelineRunResponse.class);

    // response should include the job report and pipeline output object
    assertEquals(newJobId.toString(), response.getJobReport().getId());
    assertEquals(ApiJobReport.StatusEnum.SUCCEEDED, response.getJobReport().getStatus());
    assertEquals(createdTime.toString(), response.getJobReport().getSubmitted());
    assertEquals(updatedTime.toString(), response.getJobReport().getCompleted());
    assertEquals(testOutput, response.getPipelineOutput());
    assertNull(response.getErrorReport());
  }

  @Test
  void getPipelineRunResultDoneFailed() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String jobIdString = newJobId.toString();
    String errorMessage = "test exception message";
    Integer statusCode = 500;
    PipelineRun pipelineRun = createPipelineRunCompleted(CommonPipelineRunStatusEnum.FAILED);

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
    when(pipelineRunsServiceMock.getPipelineRun(newJobId, testUser.getSubjectId()))
        .thenReturn(pipelineRun);
    when(jobServiceMock.retrieveAsyncJobResult(
            newJobId, testUser.getSubjectId(), String.class, null))
        .thenReturn(jobResult);

    MvcResult result =
        mockMvc
            .perform(
                get(String.format("/api/pipelineruns/v1/%s/result/%s", pipelineName, jobIdString)))
            .andExpect(status().isOk()) // the call itself should return a 200
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAsyncPipelineRunResponse response =
        new ObjectMapper()
            .readValue(
                result.getResponse().getContentAsString(), ApiAsyncPipelineRunResponse.class);

    // response should include the error report and no pipeline output object
    assertEquals(newJobId.toString(), response.getJobReport().getId());
    assertNull(response.getPipelineOutput());
    assertEquals(statusCode, response.getJobReport().getStatusCode());
    assertEquals(errorMessage, response.getErrorReport().getMessage());
  }

  @Test
  void getPipelineRunResultRunning() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String jobIdString = newJobId.toString();
    Integer statusCode = 202;
    PipelineRun pipelineRun = createPipelineRunRunning();
    JobApiUtils.AsyncJobResult<String> jobResult =
        new JobApiUtils.AsyncJobResult<String>()
            .jobReport(
                new ApiJobReport()
                    .id(newJobId.toString())
                    .status(ApiJobReport.StatusEnum.RUNNING)
                    .statusCode(statusCode));

    // the mocks
    when(pipelineRunsServiceMock.getPipelineRun(newJobId, testUser.getSubjectId()))
        .thenReturn(pipelineRun);
    when(jobServiceMock.retrieveAsyncJobResult(
            newJobId, testUser.getSubjectId(), String.class, null))
        .thenReturn(jobResult);

    MvcResult result =
        mockMvc
            .perform(
                get(String.format("/api/pipelineruns/v1/%s/result/%s", pipelineName, jobIdString)))
            .andExpect(status().isAccepted())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAsyncPipelineRunResponse response =
        new ObjectMapper()
            .readValue(
                result.getResponse().getContentAsString(), ApiAsyncPipelineRunResponse.class);

    // response should include the job report and no error report or pipeline output object
    assertEquals(newJobId.toString(), response.getJobReport().getId());
    assertEquals(statusCode, response.getJobReport().getStatusCode());
    assertNull(response.getPipelineOutput());
    assertNull(response.getErrorReport());
  }

  @Test
  void getPipelineRunResultNotFound() throws Exception {
    String pipelineName = PipelinesEnum.IMPUTATION_BEAGLE.getValue();
    String jobIdString = newJobId.toString();

    // the mocks
    when(pipelineRunsServiceMock.getPipelineRun(newJobId, testUser.getSubjectId()))
        .thenReturn(null);

    // the call should return a 404
    mockMvc
        .perform(get(String.format("/api/pipelineruns/v1/%s/result/%s", pipelineName, jobIdString)))
        .andExpect(status().isNotFound())
        .andExpect(
            result -> assertInstanceOf(NotFoundException.class, result.getResolvedException()));
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
        PipelineRunsApiController.getAsyncResultEndpoint(ingressConfiguration, request, jobId));
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
        PipelineRunsApiController.getAsyncResultEndpoint(ingressConfiguration, request, jobId));
  }

  // support methods

  private String createTestPipelineRunPostBody(String jobId, String description)
      throws JsonProcessingException {
    String stringifiedInputs = MockMvcUtils.convertToJsonString(testPipelineInputs);
    return String.format(
        "{\"jobControl\":{\"id\":\"%s\"},\"pipelineVersion\":\"%s\",\"pipelineInputs\":%s,\"description\":\"%s\"}",
        jobId, testPipelineVersion, stringifiedInputs, description);
  }

  /** helper method to create a PipelineRun object for a running job */
  private PipelineRun createPipelineRunRunning() {
    return new PipelineRun(
        newJobId,
        testUser.getSubjectId(),
        1L,
        TestUtils.CONTROL_WORKSPACE_ID,
        createdTime,
        updatedTime,
        CommonPipelineRunStatusEnum.RUNNING.toString(),
        TestUtils.TEST_PIPELINE_DESCRIPTION_1,
        testResultPath,
        null,
        null);
  }

  /**
   * helper method to create a PipelineRun object for a completed job. if the job status is
   * SUCCEEDED, we include an output.
   */
  private PipelineRun createPipelineRunCompleted(CommonPipelineRunStatusEnum status) {
    Boolean isSuccess = status == CommonPipelineRunStatusEnum.SUCCEEDED ? true : null;
    String output = status == CommonPipelineRunStatusEnum.SUCCEEDED ? testOutputString : null;
    return new PipelineRun(
        newJobId,
        testUser.getSubjectId(),
        1L,
        TestUtils.CONTROL_WORKSPACE_ID,
        createdTime,
        updatedTime,
        status.toString(),
        TestUtils.TEST_PIPELINE_DESCRIPTION_1,
        testResultPath,
        isSuccess,
        output);
  }
}
