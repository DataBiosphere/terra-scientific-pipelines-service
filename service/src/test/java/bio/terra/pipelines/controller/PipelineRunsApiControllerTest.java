package bio.terra.pipelines.controller;

import static bio.terra.pipelines.testutils.MockMvcUtils.getTestPipeline;
import static bio.terra.pipelines.testutils.TestUtils.buildTestResultUrl;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.*;
import static org.mockito.Mockito.doNothing;
import static org.mockito.Mockito.doThrow;
import static org.mockito.Mockito.when;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.post;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.NotFoundException;
import bio.terra.common.exception.ValidationException;
import bio.terra.common.iam.BearerTokenFactory;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.IngressConfiguration;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.app.configuration.internal.PipelinesCommonConfiguration;
import bio.terra.pipelines.app.controller.GlobalExceptionHandler;
import bio.terra.pipelines.app.controller.JobApiUtils;
import bio.terra.pipelines.app.controller.PipelineRunsApiController;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.pagination.*;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.dependencies.stairway.exception.InternalStairwayException;
import bio.terra.pipelines.generated.model.*;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
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
import java.time.Instant;
import java.time.temporal.ChronoUnit;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.http.HttpStatus;
import org.springframework.http.MediaType;
import org.springframework.http.converter.HttpMessageNotReadableException;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.bean.override.mockito.MockitoBean;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.MvcResult;
import org.springframework.test.web.servlet.result.MockMvcResultMatchers;
import org.springframework.web.bind.MethodArgumentNotValidException;

@ContextConfiguration(classes = {PipelineRunsApiController.class, GlobalExceptionHandler.class})
@WebMvcTest
class PipelineRunsApiControllerTest {
  @MockitoBean PipelinesService pipelinesServiceMock;
  @MockitoBean PipelineRunsService pipelineRunsServiceMock;
  @MockitoBean PipelineInputsOutputsService pipelineInputsOutputsServiceMock;
  @MockitoBean JobService jobServiceMock;
  @MockitoBean SamUserFactory samUserFactoryMock;
  @MockitoBean BearerTokenFactory bearerTokenFactory;
  @MockitoBean SamConfiguration samConfiguration;
  @MockitoBean IngressConfiguration ingressConfiguration;
  @MockitoBean PipelinesCommonConfiguration pipelinesCommonConfiguration;

  @Autowired private MockMvc mockMvc;

  private final SamUser testUser = MockMvcUtils.TEST_SAM_USER;
  private final int testPipelineVersion = TestUtils.TEST_PIPELINE_VERSION_1;
  private final String testPipelineToolVersion = TestUtils.TEST_TOOL_VERSION_1;
  private final Map<String, Object> testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID newJobId = TestUtils.TEST_NEW_UUID;
  private final Instant createdTime = Instant.now();
  private final Instant updatedTime = Instant.now();
  private final Map<String, String> testOutputs = TestUtils.TEST_PIPELINE_OUTPUTS;

  private final Integer testQuotaConsumed = 10;

  private final Long userDataTtlDays = 8L;

  @BeforeEach
  void beforeEach() {
    when(samConfiguration.baseUri()).thenReturn("baseSamUri");
    when(ingressConfiguration.getDomainName()).thenReturn(TestUtils.TEST_DOMAIN);
    when(samUserFactoryMock.from(any(HttpServletRequest.class), eq("baseSamUri")))
        .thenReturn(testUser);
    when(pipelinesServiceMock.getPipeline(any(PipelinesEnum.class), anyInt()))
        .thenReturn(getTestPipeline());
    when(pipelinesServiceMock.getPipelineById(anyLong())).thenReturn(getTestPipeline());
    when(pipelinesServiceMock.getPipelines()).thenReturn(List.of(getTestPipeline()));
    when(pipelinesCommonConfiguration.getUserDataTtlDays()).thenReturn(userDataTtlDays);
  }

  // preparePipelineRun tests
  // preparePipelineRun performs the following actions:
  // 1. extract user-provided information from the request
  // 2. validate the user-provided information
  // 3. call pipelineRunsService.preparePipelineRun to prepare the job
  // 4. configure a response object

  @Test
  void prepareRunImputationPipeline() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    UUID jobId = newJobId;
    String description = "description for testPrepareJobImputationPipeline";
    String postBodyAsJson =
        testPreparePipelineRunPostBody(jobId.toString(), pipelineName, description);

    Map<String, Map<String, String>> pipelineInputsWithSasUrls = new HashMap<>();
    // the contents of this doesn't matter
    testPipelineInputs.forEach(
        (key, value) -> pipelineInputsWithSasUrls.put(key, Map.of("sasUrl", value.toString())));

    // the mocks
    doNothing()
        .when(pipelineInputsOutputsServiceMock)
        .validateUserProvidedInputs(
            getTestPipeline().getPipelineInputDefinitions(), TestUtils.TEST_PIPELINE_INPUTS);
    when(pipelineRunsServiceMock.preparePipelineRun(
            getTestPipeline(),
            jobId,
            testUser.getSubjectId(),
            TestUtils.TEST_PIPELINE_INPUTS,
            description))
        .thenReturn(pipelineInputsWithSasUrls);

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post("/api/pipelineruns/v1/prepare")
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiPreparePipelineRunResponse response =
        new ObjectMapper()
            .readValue(
                result.getResponse().getContentAsString(), ApiPreparePipelineRunResponse.class);
    assertEquals(jobId, response.getJobId());
    assertEquals(pipelineInputsWithSasUrls, response.getFileInputUploadUrls());
  }

  @Test
  void preparePipelineRunMissingDescriptionOk() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    UUID jobId = newJobId;
    String stringifiedInputs = MockMvcUtils.convertToJsonString(testPipelineInputs);
    String postBodyAsJson =
        String.format(
            "{\"jobId\":\"%s\",\"pipelineName\":\"%s\",\"pipelineVersion\":\"%s\",\"pipelineInputs\":%s}",
            jobId, pipelineName, testPipelineVersion, stringifiedInputs);

    Map<String, Map<String, String>> pipelineInputsWithSasUrls = new HashMap<>();
    // the contents of this doesn't matter
    testPipelineInputs.forEach(
        (key, value) -> pipelineInputsWithSasUrls.put(key, Map.of("sasUrl", value.toString())));

    // the mocks
    doNothing()
        .when(pipelineInputsOutputsServiceMock)
        .validateUserProvidedInputs(
            getTestPipeline().getPipelineInputDefinitions(), TestUtils.TEST_PIPELINE_INPUTS);
    when(pipelineRunsServiceMock.preparePipelineRun(
            getTestPipeline(),
            jobId,
            testUser.getSubjectId(),
            TestUtils.TEST_PIPELINE_INPUTS,
            null))
        .thenReturn(pipelineInputsWithSasUrls);

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post("/api/pipelineruns/v1/prepare")
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiPreparePipelineRunResponse response =
        new ObjectMapper()
            .readValue(
                result.getResponse().getContentAsString(), ApiPreparePipelineRunResponse.class);
    assertEquals(jobId, response.getJobId());
    assertEquals(pipelineInputsWithSasUrls, response.getFileInputUploadUrls());
  }

  @Test
  void preparePipelineRunMissingMultipleRequiredFields() throws Exception {
    String stringifiedInputs = MockMvcUtils.convertToJsonString(testPipelineInputs);
    String postBodyAsJson =
        String.format(
            "{\"pipelineInputs\":%s}", // missing jobId and pipelineName and pipelineVersion
            stringifiedInputs);

    // Spring will catch the missing fields and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post("/api/pipelineruns/v1/prepare")
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
                    "Request could not be parsed or was invalid: jobId must not be null; "
                        + "pipelineName must not be null"));
  }

  @Test
  void preparePipelineRunBadPipelineInputs() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    String description = "description for testPrepareJobBadPipelineInputs";
    String postBodyAsJson =
        testPreparePipelineRunPostBody(newJobId.toString(), pipelineName, description);

    // the mocks
    doThrow(new ValidationException("some message"))
        .when(pipelineInputsOutputsServiceMock)
        .validateUserProvidedInputs(
            TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST, testPipelineInputs);

    mockMvc
        .perform(
            post("/api/pipelineruns/v1/prepare")
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isBadRequest())
        .andExpect(
            result -> assertInstanceOf(ValidationException.class, result.getResolvedException()));
  }

  // startPipelineRun tests
  // startPipelineRun performs the following actions:
  // 1. extract user-provided information from the request
  // 3. call pipelineRunsService.startPipelineRun to update pipeline_runs db and create
  // a new job in Stairway (in a transaction)
  // 4. query the pipeline_runs table for the job and configure a response object

  @Test
  void startRunImputationPipelineRunning() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    String description = "description for testCreateJobImputationPipelineRunning";
    UUID jobId = newJobId;
    String postBodyAsJson = testStartPipelineRunPostBody(jobId.toString());
    FlightMap inputParameters = StairwayTestUtils.constructCreateJobInputs(new FlightMap());
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.RUNNING, jobId, inputParameters, new FlightMap());
    PipelineRun testPipelinePrepared = getPipelineRunPreparing(description);
    PipelineRun testPipelineRun = getPipelineRunRunning();
    testPipelineRun.setDescription(description);

    ApiJobReport jobReport =
        new ApiJobReport()
            .id(flightState.getFlightId())
            .description(description)
            .status(ApiJobReport.StatusEnum.RUNNING)
            .statusCode(HttpStatus.ACCEPTED.value())
            .submitted(testPipelineRun.getCreated().toString())
            .completed(null)
            .resultURL(
                JobApiUtils.getAsyncResultEndpoint(
                    TestUtils.TEST_DOMAIN, UUID.fromString(flightState.getFlightId())));

    // the mocks
    when(pipelineRunsServiceMock.getPipelineRun(jobId, testUser.getSubjectId()))
        .thenReturn(testPipelinePrepared);
    when(pipelineRunsServiceMock.startPipelineRun(
            getTestPipeline(), jobId, testUser.getSubjectId()))
        .thenReturn(testPipelineRun);
    when(jobServiceMock.retrieveJob(jobId, testUser.getSubjectId(), PipelinesEnum.ARRAY_IMPUTATION))
        .thenReturn(flightState);
    when(jobServiceMock.retrieveAsyncJobResult(jobId, testUser.getSubjectId(), String.class, null))
        .thenReturn(new JobApiUtils.AsyncJobResult<String>().jobReport(jobReport).result(null));

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post("/api/pipelineruns/v1/start")
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isAccepted())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAsyncPipelineRunResponse response =
        new ObjectMapper()
            .readValue(
                result.getResponse().getContentAsString(), ApiAsyncPipelineRunResponse.class);
    ApiPipelineRunReport pipelineRunReportResponse = response.getPipelineRunReport();

    assertEquals(jobId.toString(), response.getJobReport().getId());
    assertEquals(buildTestResultUrl(jobId.toString()), response.getJobReport().getResultURL());
    assertEquals(ApiJobReport.StatusEnum.RUNNING, response.getJobReport().getStatus());
    assertEquals(createdTime.toString(), response.getJobReport().getSubmitted());
    assertEquals(pipelineName, pipelineRunReportResponse.getPipelineName());
    assertEquals(testPipelineVersion, pipelineRunReportResponse.getPipelineVersion());
    assertEquals(testPipelineToolVersion, pipelineRunReportResponse.getToolVersion());
    assertNull(pipelineRunReportResponse.getOutputs());
    assertNull(response.getJobReport().getCompleted());
  }

  @Test
  void startPipelineRunMissingJobControl() throws Exception {
    ApiStartPipelineRunRequestBody postBody = new ApiStartPipelineRunRequestBody();
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // Spring will catch the missing jobControl and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post("/api/pipelineruns/v1/start")
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
  void startPipelineRunMissingJobId() throws Exception {
    ApiJobControl apiJobControl = new ApiJobControl();
    ApiStartPipelineRunRequestBody postBody =
        new ApiStartPipelineRunRequestBody().jobControl(apiJobControl);
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // Spring will catch the missing job id and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post("/api/pipelineruns/v1/start")
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
  void startPipelineExpiredPipeline() throws Exception {
    String postBodyAsJson = testStartPipelineRunPostBody(newJobId.toString());
    PipelineRun expiredPreparingPipelineRun =
        getPipelineRunWithStatusAndQuotaConsumed(CommonPipelineRunStatusEnum.PREPARING, null);
    expiredPreparingPipelineRun.setCreated(
        expiredPreparingPipelineRun.getCreated().minus(userDataTtlDays + 1L, ChronoUnit.DAYS));
    when(pipelineRunsServiceMock.getPipelineRun(newJobId, testUser.getSubjectId()))
        .thenReturn(expiredPreparingPipelineRun);

    MvcResult result =
        mockMvc
            .perform(
                post("/api/pipelineruns/v1/start")
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isBadRequest())
            .andExpect(r -> assertInstanceOf(BadRequestException.class, r.getResolvedException()))
            .andReturn();

    ApiErrorReport response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiErrorReport.class);

    assertEquals(
        "Pipeline run was prepared more than %s days ago; it cannot be started"
            .formatted(userDataTtlDays),
        response.getMessage());
  }

  @Test
  void startPipelineRunBadJobId() throws Exception {
    String postBodyAsJson = testStartPipelineRunPostBody("this-is-not-a-uuid");

    // Spring will catch the non-uuid jobId and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post("/api/pipelineruns/v1/start")
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
                    "JSON parse error: Cannot deserialize value of type `java.util.UUID` "
                        + "from String \"this-is-not-a-uuid\": UUID has to be represented by "
                        + "standard 36-char representation"));
  }

  @Test
  void startPipelineRunDbError() throws Exception {
    UUID jobId = newJobId;
    String postBodyAsJson = testStartPipelineRunPostBody(jobId.toString());

    // the mocks
    when(pipelineRunsServiceMock.startPipelineRun(
            getTestPipeline(), jobId, testUser.getSubjectId()))
        .thenThrow(new RuntimeException("some message"));

    MvcResult result =
        mockMvc
            .perform(
                post("/api/pipelineruns/v1/start")
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isInternalServerError())
            .andExpect(res -> assertInstanceOf(RuntimeException.class, res.getResolvedException()))
            .andReturn();

    ApiErrorReport response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiErrorReport.class);

    assertEquals(
        "Internal server error. Please contact support if this problem persists.",
        response.getMessage());
  }

  @Test
  void startImputationRunStairwayError() throws Exception {
    UUID jobId = newJobId;
    String postBodyAsJson = testStartPipelineRunPostBody(jobId.toString());
    String description = "description for testCreateJobImputationPipelineRunning";

    // the mocks - one error that can happen is a MissingRequiredFieldException from Stairway
    when(pipelineRunsServiceMock.getPipelineRun(jobId, testUser.getSubjectId()))
        .thenReturn(getPipelineRunPreparing(description));
    when(pipelineRunsServiceMock.startPipelineRun(
            getTestPipeline(), jobId, testUser.getSubjectId()))
        .thenThrow(new InternalStairwayException("some message"));

    MvcResult result =
        mockMvc
            .perform(
                post("/api/pipelineruns/v1/start")
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(postBodyAsJson))
            .andExpect(status().isInternalServerError())
            .andExpect(
                res ->
                    assertInstanceOf(InternalStairwayException.class, res.getResolvedException()))
            .andReturn();

    ApiErrorReport response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiErrorReport.class);

    assertEquals(
        "Internal server error. Please contact support if this problem persists.",
        response.getMessage());
  }

  // getPipelineRunResult tests

  @Test
  void getPipelineRunResultDoneSuccess() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    String jobIdString = newJobId.toString();
    PipelineRun pipelineRun =
        getPipelineRunWithStatusAndQuotaConsumed(
            CommonPipelineRunStatusEnum.SUCCEEDED, testQuotaConsumed);
    ApiPipelineRunOutputs apiPipelineRunOutputs = new ApiPipelineRunOutputs();
    apiPipelineRunOutputs.putAll(testOutputs);

    // the mocks - note we don't do anything with Stairway because all our info should be in our own
    // db
    when(pipelineRunsServiceMock.getPipelineRun(newJobId, testUser.getSubjectId()))
        .thenReturn(pipelineRun);
    when(pipelineInputsOutputsServiceMock.formatPipelineRunOutputs(pipelineRun))
        .thenReturn(apiPipelineRunOutputs);

    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/pipelineruns/v1/result/%s", jobIdString)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAsyncPipelineRunResponse response =
        new ObjectMapper()
            .readValue(
                result.getResponse().getContentAsString(), ApiAsyncPipelineRunResponse.class);
    ApiPipelineRunReport pipelineRunReportResponse = response.getPipelineRunReport();

    // response should include the job report and pipeline run report including the outputs object
    assertEquals(newJobId.toString(), response.getJobReport().getId());
    assertEquals(ApiJobReport.StatusEnum.SUCCEEDED, response.getJobReport().getStatus());
    assertEquals(createdTime.toString(), response.getJobReport().getSubmitted());
    assertEquals(updatedTime.toString(), response.getJobReport().getCompleted());
    assertEquals(pipelineName, pipelineRunReportResponse.getPipelineName());
    assertEquals(testPipelineVersion, pipelineRunReportResponse.getPipelineVersion());
    assertEquals(testPipelineToolVersion, pipelineRunReportResponse.getToolVersion());
    assertEquals(testOutputs, pipelineRunReportResponse.getOutputs());
    assertEquals(
        updatedTime.plus(userDataTtlDays, ChronoUnit.DAYS).toString(),
        pipelineRunReportResponse.getOutputExpirationDate());
    assertNull(response.getErrorReport());
  }

  @Test
  void getPipelineRunResultDoneSuccessExpiredOutputs() throws Exception {
    String jobIdString = newJobId.toString();
    PipelineRun pipelineRun =
        getPipelineRunWithStatusAndQuotaConsumed(
            CommonPipelineRunStatusEnum.SUCCEEDED, testQuotaConsumed);
    // set the updated time to 1 day ago so that the outputs are expired
    pipelineRun.setUpdated(updatedTime.minus(userDataTtlDays + 1L, ChronoUnit.DAYS));
    ApiPipelineRunOutputs apiPipelineRunOutputs = new ApiPipelineRunOutputs();
    apiPipelineRunOutputs.putAll(testOutputs);

    // the mocks - note we don't do anything with Stairway because all our info should be in our db
    when(pipelineRunsServiceMock.getPipelineRun(newJobId, testUser.getSubjectId()))
        .thenReturn(pipelineRun);
    when(pipelineInputsOutputsServiceMock.formatPipelineRunOutputs(pipelineRun))
        .thenReturn(apiPipelineRunOutputs);

    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/pipelineruns/v1/result/%s", jobIdString)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAsyncPipelineRunResponse response =
        new ObjectMapper()
            .readValue(
                result.getResponse().getContentAsString(), ApiAsyncPipelineRunResponse.class);
    ApiPipelineRunReport pipelineRunReportResponse = response.getPipelineRunReport();

    // response should not include outputs because it is past the output expiration date
    assertNull(pipelineRunReportResponse.getOutputs());
  }

  @Test
  void getPipelineRunResultDoneFailed() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    String jobIdString = newJobId.toString();
    String errorMessage = "test exception message";
    Integer statusCode = 500;
    PipelineRun pipelineRun =
        getPipelineRunWithStatusAndQuotaConsumed(CommonPipelineRunStatusEnum.FAILED, null);

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
            .perform(get(String.format("/api/pipelineruns/v1/result/%s", jobIdString)))
            .andExpect(status().isOk()) // the call itself should return a 200
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAsyncPipelineRunResponse response =
        new ObjectMapper()
            .readValue(
                result.getResponse().getContentAsString(), ApiAsyncPipelineRunResponse.class);
    ApiPipelineRunReport pipelineRunReportResponse = response.getPipelineRunReport();

    // response should include the error report and pipeline run report without outputs
    assertEquals(newJobId.toString(), response.getJobReport().getId());
    assertEquals(pipelineName, pipelineRunReportResponse.getPipelineName());
    assertEquals(testPipelineVersion, pipelineRunReportResponse.getPipelineVersion());
    assertEquals(testPipelineToolVersion, pipelineRunReportResponse.getToolVersion());
    assertNull(pipelineRunReportResponse.getOutputs());
    assertEquals(statusCode, response.getJobReport().getStatusCode());
    assertEquals(errorMessage, response.getErrorReport().getMessage());
    assertNull(response.getPipelineRunReport().getOutputExpirationDate());
  }

  @Test
  void getPipelineRunResultRunning() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    String jobIdString = newJobId.toString();
    Integer statusCode = 202;
    PipelineRun pipelineRun = getPipelineRunRunning();
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
            .perform(get(String.format("/api/pipelineruns/v1/result/%s", jobIdString)))
            .andExpect(status().isAccepted())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAsyncPipelineRunResponse response =
        new ObjectMapper()
            .readValue(
                result.getResponse().getContentAsString(), ApiAsyncPipelineRunResponse.class);
    ApiPipelineRunReport pipelineRunReportResponse = response.getPipelineRunReport();

    // response should include the job report, no error report, and pipeline run report without an
    // outputs object
    assertEquals(newJobId.toString(), response.getJobReport().getId());
    assertEquals(statusCode, response.getJobReport().getStatusCode());
    assertEquals(pipelineName, pipelineRunReportResponse.getPipelineName());
    assertEquals(testPipelineVersion, pipelineRunReportResponse.getPipelineVersion());
    assertEquals(testPipelineToolVersion, pipelineRunReportResponse.getToolVersion());
    assertNull(pipelineRunReportResponse.getOutputs());
    assertNull(response.getErrorReport());
    assertNull(response.getPipelineRunReport().getOutputExpirationDate());
  }

  @Test
  void getPipelineRunResultPreparing() throws Exception {
    String jobIdString = newJobId.toString();
    String description = "description for testGetPipelineRunResultPreparing";
    PipelineRun pipelineRun = getPipelineRunPreparing(description);

    // the mocks
    when(pipelineRunsServiceMock.getPipelineRun(newJobId, testUser.getSubjectId()))
        .thenReturn(pipelineRun);

    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/pipelineruns/v1/result/%s", jobIdString)))
            .andExpect(status().isBadRequest())
            .andExpect(
                res -> assertInstanceOf(BadRequestException.class, res.getResolvedException()))
            .andReturn();

    ApiErrorReport response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiErrorReport.class);

    assertEquals(
        "Pipeline run %s is still preparing; it has to be started before you can query the result"
            .formatted(newJobId),
        response.getMessage());
  }

  @Test
  void getPipelineRunResultNotFound() throws Exception {
    String jobIdString = newJobId.toString();

    // the mocks
    when(pipelineRunsServiceMock.getPipelineRun(newJobId, testUser.getSubjectId()))
        .thenReturn(null);

    // the call should return a 404
    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/pipelineruns/v1/result/%s", jobIdString)))
            .andExpect(status().isNotFound())
            .andExpect(res -> assertInstanceOf(NotFoundException.class, res.getResolvedException()))
            .andReturn();

    ApiErrorReport response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiErrorReport.class);

    assertEquals("Pipeline run %s not found".formatted(newJobId), response.getMessage());
  }

  @Test
  void getAllPipelineRunsWithNoPageToken() throws Exception {
    int limit = 5;
    String pageToken = null;
    String preparingDescription = "preparing job";
    PipelineRun pipelineRunPreparing = getPipelineRunPreparing(preparingDescription);
    PipelineRun pipelineRunPreparingNoDescription = getPipelineRunPreparing(null);
    PipelineRun pipelineRunSucceeded =
        getPipelineRunWithStatusAndQuotaConsumed(
            CommonPipelineRunStatusEnum.SUCCEEDED, testQuotaConsumed);
    PipelineRun pipelineRunFailed =
        getPipelineRunWithStatusAndQuotaConsumed(CommonPipelineRunStatusEnum.FAILED, null);
    PageResponse<List<PipelineRun>> pageResponse =
        new PageResponse<>(
            List.of(
                pipelineRunPreparing,
                pipelineRunSucceeded,
                pipelineRunFailed,
                pipelineRunPreparingNoDescription),
            null,
            null);

    // the mocks
    when(pipelineRunsServiceMock.findPipelineRunsPaginated(
            limit, pageToken, testUser.getSubjectId()))
        .thenReturn(pageResponse);

    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/pipelineruns/v1/pipelineruns?limit=%s", limit)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiGetPipelineRunsResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiGetPipelineRunsResponse.class);

    // response should include three pipeline runs
    assertNull(response.getPageToken());
    assertEquals(4, response.getResults().size());

    // preparing run should not have a completed time
    ApiPipelineRun responsePipelineRun1 = response.getResults().get(0);
    assertEquals(pipelineRunPreparing.getStatus().name(), responsePipelineRun1.getStatus());
    assertEquals(pipelineRunPreparing.getDescription(), responsePipelineRun1.getDescription());
    assertEquals(pipelineRunPreparing.getJobId(), responsePipelineRun1.getJobId());
    assertEquals(pipelineRunPreparing.getQuotaConsumed(), responsePipelineRun1.getQuotaConsumed());
    assertEquals(getTestPipeline().getName().getValue(), responsePipelineRun1.getPipelineName());
    assertEquals(
        pipelineRunPreparing.getCreated().toString(), responsePipelineRun1.getTimeSubmitted());
    // timestamp string should be marked as UTC, i.e. end with Z
    assertTrue(responsePipelineRun1.getTimeSubmitted().endsWith("Z"));
    assertNull(responsePipelineRun1.getTimeCompleted());

    // succeeded run should have a completed time
    ApiPipelineRun responsePipelineRun2 = response.getResults().get(1);
    assertEquals(pipelineRunSucceeded.getStatus().name(), responsePipelineRun2.getStatus());
    assertEquals(pipelineRunSucceeded.getDescription(), responsePipelineRun2.getDescription());
    assertEquals(pipelineRunSucceeded.getJobId(), responsePipelineRun2.getJobId());
    assertEquals(pipelineRunSucceeded.getQuotaConsumed(), responsePipelineRun2.getQuotaConsumed());
    assertEquals(getTestPipeline().getName().getValue(), responsePipelineRun2.getPipelineName());
    assertEquals(
        pipelineRunSucceeded.getCreated().toString(), responsePipelineRun2.getTimeSubmitted());
    // timestamp string should be marked as UTC, i.e. end with Z
    assertTrue(responsePipelineRun2.getTimeSubmitted().endsWith("Z"));
    assertEquals(
        pipelineRunSucceeded.getUpdated().toString(), responsePipelineRun2.getTimeCompleted());
    // timestamp string should be marked as UTC, i.e. end with Z
    assertTrue(responsePipelineRun2.getTimeCompleted().endsWith("Z"));

    // failed run should have a completed time
    ApiPipelineRun responsePipelineRun3 = response.getResults().get(2);
    assertEquals(pipelineRunFailed.getStatus().name(), responsePipelineRun3.getStatus());
    assertEquals(pipelineRunFailed.getDescription(), responsePipelineRun3.getDescription());
    assertEquals(pipelineRunFailed.getJobId(), responsePipelineRun3.getJobId());
    assertEquals(pipelineRunFailed.getQuotaConsumed(), responsePipelineRun3.getQuotaConsumed());
    assertEquals(getTestPipeline().getName().getValue(), responsePipelineRun3.getPipelineName());
    assertEquals(
        pipelineRunFailed.getCreated().toString(), responsePipelineRun3.getTimeSubmitted());
    assertEquals(
        pipelineRunFailed.getUpdated().toString(), responsePipelineRun3.getTimeCompleted());

    // preparing run without description should not have a description
    ApiPipelineRun responsePipelineRun4 = response.getResults().get(3);
    assertNull(responsePipelineRun4.getDescription());
  }

  @Test
  void getAllPipelineRunsWithPageTokenWithNextPageOverMaxLimit() throws Exception {
    int limit = 105; // this should be limited to 100 inside of controller code
    String requestPageToken = "requestPageToken";
    String nextPageToken = "nextPageToken";
    PipelineRun pipelineRun = getPipelineRunRunning();
    PageResponse<List<PipelineRun>> pageResponse =
        new PageResponse<>(List.of(pipelineRun), null, nextPageToken);

    // mocks

    // hardcoding the limit to 100 here because the code should limit, if it didn't then the mock
    // wouldn't work and the test would fail
    when(pipelineRunsServiceMock.findPipelineRunsPaginated(
            100, requestPageToken, testUser.getSubjectId()))
        .thenReturn(pageResponse);

    MvcResult result =
        mockMvc
            .perform(
                get(
                    String.format(
                        "/api/pipelineruns/v1/pipelineruns?limit=%s&pageToken=%s",
                        limit, requestPageToken)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiGetPipelineRunsResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiGetPipelineRunsResponse.class);

    // response should one pipeline run
    assertEquals(nextPageToken, response.getPageToken());
    assertEquals(1, response.getResults().size());
    ApiPipelineRun responsePipelineRun = response.getResults().get(0);
    assertEquals(pipelineRun.getStatus().name(), responsePipelineRun.getStatus());
    assertEquals(pipelineRun.getDescription(), responsePipelineRun.getDescription());
    assertEquals(pipelineRun.getJobId(), responsePipelineRun.getJobId());
  }

  // support methods

  private String testPreparePipelineRunPostBody(
      String jobId, String pipelineName, String description) throws JsonProcessingException {
    String stringifiedInputs = MockMvcUtils.convertToJsonString(testPipelineInputs);
    return String.format(
        "{\"jobId\":\"%s\",\"pipelineName\":\"%s\",\"pipelineVersion\":\"%s\",\"pipelineInputs\":%s, \"description\":\"%s\"}",
        jobId, pipelineName, testPipelineVersion, stringifiedInputs, description);
  }

  private String testStartPipelineRunPostBody(String jobId) {
    return String.format("{\"jobControl\":{\"id\":\"%s\"}}", jobId);
  }

  /** helper method to create a PipelineRun object for a running job */
  private PipelineRun getPipelineRunPreparing(String description) {
    PipelineRun preparingPipelineRun =
        getPipelineRunWithStatusAndQuotaConsumed(CommonPipelineRunStatusEnum.PREPARING, null);
    preparingPipelineRun.setDescription(description);
    return preparingPipelineRun;
  }

  /** helper method to create a PipelineRun object for a running job */
  private PipelineRun getPipelineRunRunning() {
    return getPipelineRunWithStatusAndQuotaConsumed(CommonPipelineRunStatusEnum.RUNNING, null);
  }

  /**
   * helper method to create a PipelineRun object for a completed job, specifying the status and
   * quotaConsumed.
   */
  private PipelineRun getPipelineRunWithStatusAndQuotaConsumed(
      CommonPipelineRunStatusEnum status, Integer quotaConsumed) {
    return new PipelineRun(
        newJobId,
        testUser.getSubjectId(),
        TestUtils.TEST_PIPELINE_ID_1,
        TestUtils.TEST_TOOL_VERSION_1,
        TestUtils.CONTROL_WORKSPACE_ID,
        TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
        TestUtils.CONTROL_WORKSPACE_NAME,
        TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME,
        TestUtils.CONTROL_WORKSPACE_GOOGLE_PROJECT,
        createdTime,
        updatedTime,
        status,
        TestUtils.TEST_PIPELINE_DESCRIPTION_1,
        quotaConsumed);
  }
}
