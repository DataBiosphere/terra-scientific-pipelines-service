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
import bio.terra.pipelines.app.controller.GlobalExceptionHandler;
import bio.terra.pipelines.app.controller.JobsApiController;
import bio.terra.pipelines.db.exception.ImputationJobNotFoundException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJob;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.generated.model.ApiGetJobsResponse;
import bio.terra.pipelines.generated.model.ApiJobReport;
import bio.terra.pipelines.service.ImputationService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.testutils.MockMvcUtils;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.time.Instant;
import java.util.*;
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

@ContextConfiguration(classes = {JobsApiController.class, GlobalExceptionHandler.class})
@WebMvcTest()
class JobsApiControllerTest {
  @MockBean StairwayJobService stairwayJobServiceMock;
  @MockBean PipelinesService pipelinesServiceMock;
  @MockBean SamUserFactory samUserFactoryMock;
  @MockBean BearerTokenFactory bearerTokenFactory;
  @MockBean SamConfiguration samConfiguration;
  @MockBean SamService samService;
  @MockBean ImputationService imputationService;

  @Autowired private MockMvc mockMvc;
  private final SamUser testUser = MockMvcUtils.TEST_SAM_USER;
  private final String testUserId = testUser.getSubjectId();
  private final String pipelineId = TestUtils.TEST_PIPELINE_ID_1;
  private final String pipelineVersion = TestUtils.TEST_PIPELINE_VERSION_1;
  // should be updated once we do more thinking on what this will look like
  private final Object pipelineInputs = Collections.emptyMap();
  private final Instant timeSubmittedOne = Instant.now();
  private final Instant timeSubmittedTwo = Instant.now();
  private final Instant timeCompletedOne = Instant.now();
  private final Instant timeCompletedTwo = Instant.now();
  private final UUID jobIdOkDone = TestUtils.TEST_NEW_UUID;
  private final UUID secondJobId = UUID.randomUUID();
  private final FlightMap createJobInputParameters =
      StairwayTestUtils.constructCreateJobInputs(
          pipelineId, pipelineVersion, testUserId, pipelineInputs);
  private final FlightMap createJobWorkingMap = new FlightMap();
  private final FlightState flightStateDoneSuccess =
      StairwayTestUtils.constructFlightStateWithStatusAndId(
          FlightStatus.SUCCESS,
          jobIdOkDone,
          createJobInputParameters,
          createJobWorkingMap,
          timeSubmittedOne,
          timeCompletedOne);
  private final FlightState secondFlightStateDoneSuccess =
      StairwayTestUtils.constructFlightStateWithStatusAndId(
          FlightStatus.SUCCESS,
          secondJobId,
          createJobInputParameters,
          createJobWorkingMap,
          timeSubmittedTwo,
          timeCompletedTwo);

  private final EnumeratedJob jobDoneSuccess =
      new EnumeratedJob().flightState(flightStateDoneSuccess);
  private final EnumeratedJob secondJobDoneSuccess =
      new EnumeratedJob().flightState(secondFlightStateDoneSuccess);

  @BeforeEach
  void beforeEach() {
    when(samUserFactoryMock.from(any(HttpServletRequest.class), any())).thenReturn(testUser);
  }

  @Test
  void testGetJobOk() throws Exception {
    when(stairwayJobServiceMock.retrieveJob(jobIdOkDone, testUserId))
        .thenReturn(flightStateDoneSuccess);

    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/job/v1alpha1/jobs/%s", jobIdOkDone)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiJobReport response =
        new ObjectMapper().readValue(result.getResponse().getContentAsString(), ApiJobReport.class);

    // you could compare other fields here too beyond the id, if wanted
    assertEquals(jobIdOkDone.toString(), response.getId());
  }

  @Test
  void testGetErrorJobOk() throws Exception {
    UUID jobId = UUID.randomUUID();
    FlightState flightStateDoneError =
        StairwayTestUtils.constructFlightStateWithStatusAndId(
            FlightStatus.ERROR,
            jobId,
            createJobInputParameters,
            createJobWorkingMap,
            timeSubmittedOne,
            timeCompletedOne);

    when(stairwayJobServiceMock.retrieveJob(jobId, testUserId)).thenReturn(flightStateDoneError);

    // even though the job itself failed, it completed successfully so the status code should be 200
    // (ok)
    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/job/v1alpha1/jobs/%s", jobId)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiJobReport response =
        new ObjectMapper().readValue(result.getResponse().getContentAsString(), ApiJobReport.class);

    // you could compare other fields here too beyond the id, if wanted
    assertEquals(jobId.toString(), response.getId());
  }

  @Test
  void testGetJobNotFound() throws Exception {
    UUID badJobId = UUID.randomUUID();
    when(stairwayJobServiceMock.retrieveJob(badJobId, testUserId))
        .thenThrow(new ImputationJobNotFoundException("some message"));

    mockMvc
        .perform(get(String.format("/api/job/v1alpha1/jobs/%s", badJobId)))
        .andExpect(status().isNotFound());
  }

  @Test
  void testGetMultipleJobs() throws Exception {
    EnumeratedJobs bothJobs =
        new EnumeratedJobs().results(List.of(jobDoneSuccess, secondJobDoneSuccess)).totalResults(2);

    // the mocks
    when(pipelinesServiceMock.pipelineExists(pipelineId)).thenReturn(true);
    when(stairwayJobServiceMock.enumerateJobs(testUser.getSubjectId(), 10, null, null))
        .thenReturn(bothJobs);

    MvcResult result =
        mockMvc
            .perform(get("/api/job/v1alpha1/jobs"))
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
