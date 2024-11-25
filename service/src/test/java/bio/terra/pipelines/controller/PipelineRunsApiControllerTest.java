package bio.terra.pipelines.controller;

import static bio.terra.pipelines.testutils.MockMvcUtils.getTestPipeline;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.ArgumentMatchers.*;
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
import java.time.LocalDateTime;
import java.util.HashMap;
import java.util.List;
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
  @MockBean PipelineInputsOutputsService pipelineInputsOutputsServiceMock;
  @MockBean JobService jobServiceMock;
  @MockBean SamUserFactory samUserFactoryMock;
  @MockBean BearerTokenFactory bearerTokenFactory;
  @MockBean SamConfiguration samConfiguration;
  @MockBean IngressConfiguration ingressConfiguration;

  @Autowired private MockMvc mockMvc;

  private final SamUser testUser = MockMvcUtils.TEST_SAM_USER;
  private final int testPipelineVersion = TestUtils.TEST_PIPELINE_VERSION_1;
  private final String testPipelineWdlMethodVersion = TestUtils.TEST_WDL_METHOD_VERSION_1;
  private final Map<String, Object> testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID newJobId = TestUtils.TEST_NEW_UUID;
  private final LocalDateTime createdTime = LocalDateTime.now();
  private final LocalDateTime updatedTime = LocalDateTime.now();
  private final String testResultPath = TestUtils.TEST_RESULT_URL;
  private final Map<String, String> testOutputs = TestUtils.TEST_PIPELINE_OUTPUTS;

  @BeforeEach
  void beforeEach() {
    when(samConfiguration.baseUri()).thenReturn("baseSamUri");
    when(ingressConfiguration.getDomainName()).thenReturn(TestUtils.TEST_DOMAIN);
    when(samUserFactoryMock.from(any(HttpServletRequest.class), eq("baseSamUri")))
        .thenReturn(testUser);
    when(pipelinesServiceMock.getPipeline(any(PipelinesEnum.class))).thenReturn(getTestPipeline());
    when(pipelinesServiceMock.getPipelineById(anyLong())).thenReturn(getTestPipeline());
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
    String postBodyAsJson = testPreparePipelineRunPostBody(jobId.toString());

    Map<String, Map<String, String>> pipelineInputsWithSasUrls = new HashMap<>();
    // the contents of this doesn't matter
    testPipelineInputs.forEach(
        (key, value) -> pipelineInputsWithSasUrls.put(key, Map.of("sasUrl", value.toString())));

    // the mocks
    doNothing()
        .when(pipelinesServiceMock)
        .validateUserProvidedInputs(
            getTestPipeline().getPipelineInputDefinitions(), TestUtils.TEST_PIPELINE_INPUTS);
    when(pipelineRunsServiceMock.preparePipelineRun(
            getTestPipeline(), jobId, testUser.getSubjectId(), TestUtils.TEST_PIPELINE_INPUTS))
        .thenReturn(pipelineInputsWithSasUrls);

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelineruns/v1/%s/prepare", pipelineName))
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
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    String stringifiedInputs = MockMvcUtils.convertToJsonString(testPipelineInputs);
    String postBodyAsJson =
        String.format(
            "{\"pipelineInputs\":%s}", // missing jobId and pipelineVersion
            stringifiedInputs);

    // Spring will catch the missing fields and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelineruns/v1/%s/prepare", pipelineName))
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
                        + "pipelineVersion must not be null"));
  }

  @Test
  void preparePipelineRunBadPipelineInputs() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    String postBodyAsJson = testPreparePipelineRunPostBody(newJobId.toString());

    // the mocks
    doThrow(new ValidationException("some message"))
        .when(pipelinesServiceMock)
        .validateUserProvidedInputs(
            TestUtils.TEST_PIPELINE_INPUTS_DEFINITION_LIST, testPipelineInputs);

    mockMvc
        .perform(
            post(String.format("/api/pipelineruns/v1/%s/prepare", pipelineName))
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
    String postBodyAsJson = testStartPipelineRunPostBody(jobId.toString(), description);
    FlightMap inputParameters = StairwayTestUtils.constructCreateJobInputs(new FlightMap());
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.RUNNING, jobId, inputParameters, new FlightMap());
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
            .resultURL(testPipelineRun.getResultUrl());

    // the mocks
    when(pipelineRunsServiceMock.startPipelineRun(
            getTestPipeline(),
            jobId,
            testUser.getSubjectId(),
            description,
            "https://some-teaspoons-domain.com/result/deadbeef-dead-beef-aaaa-beefdeadbeef"))
        .thenReturn(testPipelineRun);
    when(jobServiceMock.retrieveJob(jobId, testUser.getSubjectId(), PipelinesEnum.ARRAY_IMPUTATION))
        .thenReturn(flightState);
    when(jobServiceMock.retrieveAsyncJobResult(jobId, testUser.getSubjectId(), String.class, null))
        .thenReturn(new JobApiUtils.AsyncJobResult<String>().jobReport(jobReport).result(null));

    // make the call
    MvcResult result =
        mockMvc
            .perform(
                post(String.format("/api/pipelineruns/v1/%s/start", pipelineName))
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
    assertEquals(testResultPath, response.getJobReport().getResultURL());
    assertEquals(ApiJobReport.StatusEnum.RUNNING, response.getJobReport().getStatus());
    assertEquals(createdTime.toString(), response.getJobReport().getSubmitted());
    assertEquals(pipelineName, pipelineRunReportResponse.getPipelineName());
    assertEquals(testPipelineVersion, pipelineRunReportResponse.getPipelineVersion());
    assertEquals(testPipelineWdlMethodVersion, pipelineRunReportResponse.getWdlMethodVersion());
    assertNull(pipelineRunReportResponse.getOutputs());
    assertNull(response.getJobReport().getCompleted());
  }

  @Test
  void startPipelineRunMissingJobControl() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    ApiStartPipelineRunRequestBody postBody =
        new ApiStartPipelineRunRequestBody()
            .description("description for testCreateJobMissingJobId");
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // Spring will catch the missing jobControl and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelineruns/v1/%s/start", pipelineName))
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
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    ApiJobControl apiJobControl = new ApiJobControl();
    ApiStartPipelineRunRequestBody postBody =
        new ApiStartPipelineRunRequestBody()
            .jobControl(apiJobControl)
            .description("description for testCreateJobMissingJobId");
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // Spring will catch the missing job id and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelineruns/v1/%s/start", pipelineName))
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
  void startPipelineRunBadJobId() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    String postBodyAsJson =
        testStartPipelineRunPostBody(
            "this-is-not-a-uuid", "description for testCreateJobMissingJobId");

    // Spring will catch the non-uuid jobId and invoke the GlobalExceptionHandler
    // before it gets to the controller
    mockMvc
        .perform(
            post(String.format("/api/pipelineruns/v1/%s/start", pipelineName))
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
  void startPipelineRunDbError() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    String description = "description for startPipelineRunDbError";
    UUID jobId = newJobId;
    String postBodyAsJson = testStartPipelineRunPostBody(jobId.toString(), description);

    // the mocks
    when(pipelineRunsServiceMock.startPipelineRun(
            getTestPipeline(), jobId, testUser.getSubjectId(), description, testResultPath))
        .thenThrow(new RuntimeException("some message"));

    mockMvc
        .perform(
            post(String.format("/api/pipelineruns/v1/%s/start", pipelineName))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isInternalServerError())
        .andExpect(
            result -> assertInstanceOf(RuntimeException.class, result.getResolvedException()));
  }

  @Test
  void startImputationRunStairwayError() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    String description = "description for startImputationJobStairwayError";
    UUID jobId = newJobId;
    String postBodyAsJson = testStartPipelineRunPostBody(jobId.toString(), description);
    String resultPath = "https://" + TestUtils.TEST_DOMAIN + "/result/" + jobId;

    // the mocks - one error that can happen is a MissingRequiredFieldException from Stairway
    when(pipelineRunsServiceMock.startPipelineRun(
            getTestPipeline(), jobId, testUser.getSubjectId(), description, resultPath))
        .thenThrow(new InternalStairwayException("some message"));

    mockMvc
        .perform(
            post(String.format("/api/pipelineruns/v1/%s/start", pipelineName))
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
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    String jobIdString = newJobId.toString();
    PipelineRun pipelineRun = getPipelineRunWithStatus(CommonPipelineRunStatusEnum.SUCCEEDED);
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
            .perform(
                get(String.format("/api/pipelineruns/v1/%s/result/%s", pipelineName, jobIdString)))
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
    assertEquals(testPipelineWdlMethodVersion, pipelineRunReportResponse.getWdlMethodVersion());
    assertEquals(testOutputs, pipelineRunReportResponse.getOutputs());
    assertNull(response.getErrorReport());
  }

  @Test
  void getPipelineRunResultDoneFailed() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
    String jobIdString = newJobId.toString();
    String errorMessage = "test exception message";
    Integer statusCode = 500;
    PipelineRun pipelineRun = getPipelineRunWithStatus(CommonPipelineRunStatusEnum.FAILED);

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
    ApiPipelineRunReport pipelineRunReportResponse = response.getPipelineRunReport();

    // response should include the error report and pipeline run report without outputs
    assertEquals(newJobId.toString(), response.getJobReport().getId());
    assertEquals(pipelineName, pipelineRunReportResponse.getPipelineName());
    assertEquals(testPipelineVersion, pipelineRunReportResponse.getPipelineVersion());
    assertEquals(testPipelineWdlMethodVersion, pipelineRunReportResponse.getWdlMethodVersion());
    assertNull(pipelineRunReportResponse.getOutputs());
    assertEquals(statusCode, response.getJobReport().getStatusCode());
    assertEquals(errorMessage, response.getErrorReport().getMessage());
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
            .perform(
                get(String.format("/api/pipelineruns/v1/%s/result/%s", pipelineName, jobIdString)))
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
    assertEquals(testPipelineWdlMethodVersion, pipelineRunReportResponse.getWdlMethodVersion());
    assertNull(pipelineRunReportResponse.getOutputs());
    assertNull(response.getErrorReport());
  }

  @Test
  void getPipelineRunResultNotFound() throws Exception {
    String pipelineName = PipelinesEnum.ARRAY_IMPUTATION.getValue();
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

  @Test
  void getAllPipelineRunsWithNoPageToken() throws Exception {
    int limit = 5;
    String pageToken = null;
    PipelineRun pipelineRunPreparing = getPipelineRunPreparing();
    PipelineRun pipelineRunSucceeded =
        getPipelineRunWithStatus(CommonPipelineRunStatusEnum.SUCCEEDED);
    PipelineRun pipelineRunFailed = getPipelineRunWithStatus(CommonPipelineRunStatusEnum.FAILED);
    PageResponse<List<PipelineRun>> pageResponse =
        new PageResponse<>(
            List.of(pipelineRunPreparing, pipelineRunSucceeded, pipelineRunFailed), null, null);

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
    assertEquals(3, response.getResults().size());

    // preparing run should not have a completed time
    ApiPipelineRun responsePipelineRun1 = response.getResults().get(0);
    assertEquals(pipelineRunPreparing.getStatus().name(), responsePipelineRun1.getStatus());
    assertEquals(pipelineRunPreparing.getDescription(), responsePipelineRun1.getDescription());
    assertEquals(pipelineRunPreparing.getJobId(), responsePipelineRun1.getJobId());
    assertEquals(
        pipelineRunPreparing.getCreated().toString(), responsePipelineRun1.getTimeSubmitted());
    assertNull(responsePipelineRun1.getTimeCompleted());

    // succeeded run should have a completed time
    ApiPipelineRun responsePipelineRun2 = response.getResults().get(1);
    assertEquals(pipelineRunSucceeded.getStatus().name(), responsePipelineRun2.getStatus());
    assertEquals(pipelineRunSucceeded.getDescription(), responsePipelineRun2.getDescription());
    assertEquals(pipelineRunSucceeded.getJobId(), responsePipelineRun2.getJobId());
    assertEquals(
        pipelineRunSucceeded.getCreated().toString(), responsePipelineRun2.getTimeSubmitted());
    assertEquals(
        pipelineRunSucceeded.getUpdated().toString(), responsePipelineRun2.getTimeCompleted());

    // failed run should have a completed time
    ApiPipelineRun responsePipelineRun3 = response.getResults().get(2);
    assertEquals(pipelineRunFailed.getStatus().name(), responsePipelineRun3.getStatus());
    assertEquals(pipelineRunFailed.getDescription(), responsePipelineRun3.getDescription());
    assertEquals(pipelineRunFailed.getJobId(), responsePipelineRun3.getJobId());
    assertEquals(
        pipelineRunFailed.getCreated().toString(), responsePipelineRun3.getTimeSubmitted());
    assertEquals(
        pipelineRunFailed.getUpdated().toString(), responsePipelineRun3.getTimeCompleted());
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

  private String testPreparePipelineRunPostBody(String jobId) throws JsonProcessingException {
    String stringifiedInputs = MockMvcUtils.convertToJsonString(testPipelineInputs);
    return String.format(
        "{\"jobId\":\"%s\",\"pipelineVersion\":\"%s\",\"pipelineInputs\":%s}",
        jobId, testPipelineVersion, stringifiedInputs);
  }

  private String testStartPipelineRunPostBody(String jobId, String description) {
    return String.format(
        "{\"jobControl\":{\"id\":\"%s\"},\"description\":\"%s\"}", jobId, description);
  }

  /** helper method to create a PipelineRun object for a running job */
  private PipelineRun getPipelineRunPreparing() {
    return new PipelineRun(
        newJobId,
        testUser.getSubjectId(),
        TestUtils.TEST_PIPELINE_ID_1,
        TestUtils.TEST_WDL_METHOD_VERSION_1,
        TestUtils.CONTROL_WORKSPACE_ID,
        TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
        TestUtils.CONTROL_WORKSPACE_NAME,
        TestUtils.CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME,
        TestUtils.CONTROL_WORKSPACE_GOOGLE_PROJECT,
        createdTime,
        updatedTime,
        CommonPipelineRunStatusEnum.PREPARING,
        TestUtils.TEST_PIPELINE_DESCRIPTION_1,
        testResultPath);
  }

  /** helper method to create a PipelineRun object for a running job */
  private PipelineRun getPipelineRunRunning() {
    return new PipelineRun(
        newJobId,
        testUser.getSubjectId(),
        TestUtils.TEST_PIPELINE_ID_1,
        TestUtils.TEST_WDL_METHOD_VERSION_1,
        TestUtils.CONTROL_WORKSPACE_ID,
        TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
        TestUtils.CONTROL_WORKSPACE_NAME,
        TestUtils.CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME,
        TestUtils.CONTROL_WORKSPACE_GOOGLE_PROJECT,
        createdTime,
        updatedTime,
        CommonPipelineRunStatusEnum.RUNNING,
        TestUtils.TEST_PIPELINE_DESCRIPTION_1,
        testResultPath);
  }

  /** helper method to create a PipelineRun object for a completed job. */
  private PipelineRun getPipelineRunWithStatus(CommonPipelineRunStatusEnum status) {
    return new PipelineRun(
        newJobId,
        testUser.getSubjectId(),
        TestUtils.TEST_PIPELINE_ID_1,
        TestUtils.TEST_WDL_METHOD_VERSION_1,
        TestUtils.CONTROL_WORKSPACE_ID,
        TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT,
        TestUtils.CONTROL_WORKSPACE_NAME,
        TestUtils.CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME,
        TestUtils.CONTROL_WORKSPACE_GOOGLE_PROJECT,
        createdTime,
        updatedTime,
        status,
        TestUtils.TEST_PIPELINE_DESCRIPTION_1,
        testResultPath);
  }
}
