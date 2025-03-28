package bio.terra.pipelines.controller;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.when;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.common.iam.BearerTokenFactory;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.app.controller.GlobalExceptionHandler;
import bio.terra.pipelines.app.controller.JobsApiController;
import bio.terra.pipelines.db.exception.ImputationJobNotFoundException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.dependencies.stairway.exception.JobUnauthorizedException;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.generated.model.ApiGetJobsResponse;
import bio.terra.pipelines.generated.model.ApiJobReport;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.testutils.MockMvcUtils;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import com.fasterxml.jackson.databind.ObjectMapper;
import jakarta.servlet.http.HttpServletRequest;
import java.util.*;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.http.MediaType;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.bean.override.mockito.MockitoBean;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.MvcResult;

@ContextConfiguration(classes = {JobsApiController.class, GlobalExceptionHandler.class})
@WebMvcTest()
class JobsApiControllerTest {
  @MockitoBean JobService jobServiceMock;
  @MockitoBean PipelinesService pipelinesServiceMock;
  @MockitoBean SamUserFactory samUserFactoryMock;
  @MockitoBean BearerTokenFactory bearerTokenFactory;
  @MockitoBean SamConfiguration samConfiguration;
  @MockitoBean SamService samService;

  @Autowired private MockMvc mockMvc;
  private final SamUser testUser = MockMvcUtils.TEST_SAM_USER;
  private final String testUserId = testUser.getSubjectId();

  @BeforeEach
  void beforeEach() {
    when(samConfiguration.baseUri()).thenReturn("baseSamUri");
    when(samUserFactoryMock.from(any(HttpServletRequest.class), eq("baseSamUri")))
        .thenReturn(testUser);
  }

  @Test
  void getJobOk() throws Exception {
    UUID jobId = TestUtils.TEST_NEW_UUID;
    FlightState flightState = StairwayTestUtils.FLIGHT_STATE_DONE_SUCCESS_1;

    when(jobServiceMock.retrieveJob(jobId, testUserId)).thenReturn(flightState);

    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/job/v1/jobs/%s", jobId)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiJobReport response =
        new ObjectMapper().readValue(result.getResponse().getContentAsString(), ApiJobReport.class);

    // you could compare other fields here too beyond the id, if wanted
    assertEquals(jobId.toString(), response.getId());
  }

  @Test
  void getErrorJobOk() throws Exception {
    UUID jobId = TestUtils.TEST_NEW_UUID;
    FlightState flightState =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.ERROR,
            jobId,
            StairwayTestUtils.CREATE_JOB_INPUT_PARAMS,
            StairwayTestUtils.EMPTY_WORKING_MAP,
            StairwayTestUtils.TIME_SUBMITTED_1,
            StairwayTestUtils.TIME_COMPLETED_1);
    flightState.setException(new Exception("Test exception"));

    when(jobServiceMock.retrieveJob(jobId, testUserId)).thenReturn(flightState);

    // even though the job itself failed, it completed successfully so the status code should be 200
    // (ok)
    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/job/v1/jobs/%s", jobId)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiJobReport response =
        new ObjectMapper().readValue(result.getResponse().getContentAsString(), ApiJobReport.class);

    // you could compare other fields here too beyond the id, if wanted
    assertEquals(jobId.toString(), response.getId());
  }

  @Test
  void getJobNotFound() throws Exception {
    UUID badJobId = UUID.randomUUID();
    when(jobServiceMock.retrieveJob(badJobId, testUserId))
        .thenThrow(new ImputationJobNotFoundException("some message"));

    mockMvc
        .perform(get(String.format("/api/job/v1/jobs/%s", badJobId)))
        .andExpect(status().isNotFound());
  }

  @Test
  void getJobNoAccess() throws Exception {
    UUID badJobId = UUID.randomUUID();
    when(jobServiceMock.retrieveJob(badJobId, testUserId))
        .thenThrow(new JobUnauthorizedException("some message"));

    mockMvc
        .perform(get(String.format("/api/job/v1/jobs/%s", badJobId)))
        .andExpect(status().isForbidden());
  }

  @Test
  void getJobBadId() throws Exception {
    String badJobId = "not-a-uuid";

    mockMvc
        .perform(get(String.format("/api/job/v1/jobs/%s", badJobId)))
        .andExpect(status().isBadRequest());
  }

  @Test
  void getMultipleJobs() throws Exception {
    EnumeratedJobs bothJobs = StairwayTestUtils.ENUMERATED_JOBS;

    // the mocks
    when(jobServiceMock.enumerateJobs(testUser.getSubjectId(), 10, null, null))
        .thenReturn(bothJobs);

    MvcResult result =
        mockMvc
            .perform(get("/api/job/v1/jobs"))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    // Now that we have the result object, we should further validate the contents of the string
    // by reconstituting the response object from the json
    ApiGetJobsResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiGetJobsResponse.class);

    // should be the same number of items as what jobsServiceMock returns
    assertEquals(bothJobs.getTotalResults(), response.getTotalResults());

    // The ids should all match what was returned from jobsServiceMock
    for (int i = 0; i < response.getTotalResults(); ++i) {
      String rawJobId = bothJobs.getResults().get(i).getFlightState().getFlightId();
      String responseJobId = response.getResults().get(i).getId();
      assertEquals(rawJobId, responseJobId);
    }
  }

  @Test
  void getAllJobsOverMaxLimit() throws Exception {
    int limit = 105; // this should be limited to 100 inside of controller code
    EnumeratedJobs bothJobs = StairwayTestUtils.ENUMERATED_JOBS;

    // mocks

    // hardcoding the limit to 100 here because the code should limit, if it didn't then the mock
    // wouldn't work and the test would fail
    when(jobServiceMock.enumerateJobs(testUser.getSubjectId(), 100, null, null))
        .thenReturn(bothJobs);

    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/job/v1/jobs?limit=%s", limit)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    // Now that we have the result object, we should further validate the contents of the string
    // by reconstituting the response object from the json
    ApiGetJobsResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiGetJobsResponse.class);

    // should be the same number of items as what jobsServiceMock returns
    assertEquals(bothJobs.getTotalResults(), response.getTotalResults());

    // The ids should all match what was returned from jobsServiceMock
    for (int i = 0; i < response.getTotalResults(); ++i) {
      String rawJobId = bothJobs.getResults().get(i).getFlightState().getFlightId();
      String responseJobId = response.getResults().get(i).getId();
      assertEquals(rawJobId, responseJobId);
    }
  }
}
