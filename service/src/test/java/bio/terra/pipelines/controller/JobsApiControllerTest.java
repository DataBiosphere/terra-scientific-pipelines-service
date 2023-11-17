package bio.terra.pipelines.controller;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.post;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.common.exception.ApiException;
import bio.terra.common.iam.BearerTokenFactory;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.app.controller.GlobalExceptionHandler;
import bio.terra.pipelines.app.controller.JobsApiController;
<<<<<<< HEAD
=======
import bio.terra.pipelines.configuration.external.SamConfiguration;
>>>>>>> 655b57e (folder reconfiguration)
import bio.terra.pipelines.db.entities.Job;
import bio.terra.pipelines.db.exception.JobNotFoundException;
import bio.terra.pipelines.db.exception.PipelineNotFoundException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.generated.model.ApiGetJobResponse;
import bio.terra.pipelines.generated.model.ApiGetJobsResponse;
import bio.terra.pipelines.generated.model.ApiPostJobRequestBody;
import bio.terra.pipelines.service.JobsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.testutils.MockMvcUtils;
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
  @MockBean JobsService jobsServiceMock;
  @MockBean PipelinesService pipelinesServiceMock;
  @MockBean SamUserFactory samUserFactoryMock;
  @MockBean BearerTokenFactory bearerTokenFactory;
  @MockBean SamConfiguration samConfiguration;
  @MockBean SamService samService;

  @Autowired private MockMvc mockMvc;
  private final SamUser testUser = MockMvcUtils.TEST_SAM_USER;

  private final String pipelineId = "imputation";
  private final String pipelineVersion = "TestVersion";
  // should be updated once we do more thinking on what this will look like
  private final Object pipelineInputs = Collections.emptyMap();
  private final Instant timestamp = Instant.now();
  private final UUID jobIdOkDone = UUID.randomUUID();
  private final UUID secondJobId = UUID.randomUUID();
  private final Job jobOkDone =
      new Job(
          jobIdOkDone,
          testUser.getSubjectId(),
          pipelineId,
          "v0",
          timestamp,
          timestamp,
          "COMPLETED");
  private final Job secondJob =
      new Job(
          secondJobId,
          testUser.getSubjectId(),
          pipelineId,
          "v0",
          timestamp,
          timestamp,
          "COMPLETED");

  @BeforeEach
  void beforeEach() {
    when(samUserFactoryMock.from(any(HttpServletRequest.class), any())).thenReturn(testUser);
  }

  @Test
  void testGetJobOk() throws Exception {

    when(jobsServiceMock.getJob(testUser.getSubjectId(), pipelineId, jobIdOkDone))
        .thenReturn(jobOkDone);

    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/jobs/v1alpha1/%s/%s", pipelineId, jobIdOkDone)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiGetJobResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiGetJobResponse.class);
    // you could compare other fields here too beyond the id, if wanted
    assertEquals(jobIdOkDone.toString(), response.getJobId());
  }

  @Test
  void testGetJobNotFound() throws Exception {
    UUID badJobId = UUID.randomUUID();
    when(jobsServiceMock.getJob(testUser.getSubjectId(), pipelineId, badJobId))
        .thenThrow(new JobNotFoundException("some message"));

    mockMvc
        .perform(get(String.format("/api/jobs/v1alpha1/%s/%s", pipelineId, badJobId)))
        .andExpect(status().isNotFound());
  }

  @Test
  void testCreateJobGoodPipeline() throws Exception {
    // This makes the body of the post... which is a lot for very little
    ApiPostJobRequestBody postBody =
        new ApiPostJobRequestBody().pipelineVersion(pipelineVersion).pipelineInputs(pipelineInputs);
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    UUID fakeJobId = UUID.randomUUID();

    // the mocks
    when(jobsServiceMock.createJob(
            testUser.getSubjectId(), pipelineId, pipelineVersion, pipelineInputs))
        .thenReturn(fakeJobId);
    when(pipelinesServiceMock.pipelineExists(pipelineId)).thenReturn(true);

    // the crafting the expected response json
    Map<String, UUID> expectedResponseMap = new HashMap<>();
    expectedResponseMap.put("jobId", fakeJobId);
    String expectedResponseJson = MockMvcUtils.convertToJsonString(expectedResponseMap);
    mockMvc
        .perform(
            post(String.format("/api/jobs/v1alpha1/%s", pipelineId))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isOk())
        .andExpect(content().contentType(MediaType.APPLICATION_JSON))
        .andExpect(content().string(expectedResponseJson));
  }

  @Test
  void testCreateJobBadPipeline() throws Exception {
    // This makes the body of the post... which is a lot for very little
    ApiPostJobRequestBody postBody =
        new ApiPostJobRequestBody().pipelineVersion(pipelineVersion).pipelineInputs(pipelineInputs);
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // the mocks
    when(pipelinesServiceMock.pipelineExists(pipelineId)).thenReturn(false);

    mockMvc
        .perform(
            post(String.format("/api/jobs/v1alpha1/%s", pipelineId))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isNotFound())
        .andExpect(
            result ->
                assertTrue(result.getResolvedException() instanceof PipelineNotFoundException));
  }

  @Test
  void testCreateJobWriteFail() throws Exception {
    // This makes the body of the post... which is a lot for very little
    ApiPostJobRequestBody postBody =
        new ApiPostJobRequestBody().pipelineVersion(pipelineVersion).pipelineInputs(pipelineInputs);
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    // the mocks - if createJob repeatedly fails to write to the database, it returns null
    when(pipelinesServiceMock.pipelineExists(pipelineId)).thenReturn(true);
    when(jobsServiceMock.createJob(
            testUser.getSubjectId(), pipelineId, pipelineVersion, pipelineInputs))
        .thenReturn(null);

    mockMvc
        .perform(
            post(String.format("/api/jobs/v1alpha1/%s", pipelineId))
                .contentType(MediaType.APPLICATION_JSON)
                .content(postBodyAsJson))
        .andExpect(status().isInternalServerError())
        .andExpect(result -> assertTrue(result.getResolvedException() instanceof ApiException));
  }

  @Test
  void testGetMultipleJobs() throws Exception {
    List<Job> bothJobs = List.of(jobOkDone, secondJob);

    // the mocks
    when(pipelinesServiceMock.pipelineExists(pipelineId)).thenReturn(true);
    when(jobsServiceMock.getJobs(testUser.getSubjectId(), pipelineId)).thenReturn(bothJobs);

    MvcResult result =
        mockMvc
            .perform(get(String.format("/api/jobs/v1alpha1/%s", pipelineId)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    // Now that we have the result object, we should further validate the contents of the string
    // by reconstituting the response object from the json
    ApiGetJobsResponse response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiGetJobsResponse.class);

    // should be the same number of items as what jobsServiceMock returns
    assertEquals(bothJobs.size(), response.size());

    // The ids should all match what was returned from jobsServiceMock
    for (int i = 0; i < response.size(); ++i) {
      String rawJobId = bothJobs.get(i).getJobId().toString();
      String responseJobId = response.get(i).getJobId();
      assertEquals(rawJobId, responseJobId);
    }
  }
}
