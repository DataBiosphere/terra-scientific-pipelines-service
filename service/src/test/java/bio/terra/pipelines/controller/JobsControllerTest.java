package bio.terra.pipelines.controller;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.post;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.common.iam.BearerToken;
import bio.terra.common.iam.BearerTokenFactory;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.controller.GlobalExceptionHandler;
import bio.terra.pipelines.app.controller.JobsApiController;
import bio.terra.pipelines.config.SamConfiguration;
import bio.terra.pipelines.db.exception.JobNotFoundException;
import bio.terra.pipelines.db.exception.PipelineNotFoundException;
import bio.terra.pipelines.generated.model.ApiGetJobsResponse;
import bio.terra.pipelines.generated.model.ApiPostJobRequestBody;
import bio.terra.pipelines.iam.SamService;
import bio.terra.pipelines.service.JobsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.service.model.Job;
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
class JobsControllerTest {
  @MockBean JobsService jobsServiceMock;
  @MockBean PipelinesService pipelinesServiceMock;
  @MockBean SamUserFactory samUserFactoryMock;
  @MockBean BearerTokenFactory bearerTokenFactory;
  @MockBean SamConfiguration samConfiguration;
  @MockBean SamService samService;

  @Autowired private MockMvc mockMvc;
  private final SamUser testUser =
      new SamUser(
          "test@email",
          UUID.randomUUID().toString(),
          new BearerToken(UUID.randomUUID().toString()));

  private final String pipelineId = "imputation";
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
          Optional.of(timestamp),
          "COMPLETED");
  private final Job secondJob =
      new Job(
          secondJobId,
          testUser.getSubjectId(),
          pipelineId,
          "v0",
          timestamp,
          Optional.of(timestamp),
          "COMPLETED");

  @BeforeEach
  void beforeEach() {
    when(samUserFactoryMock.from(any(HttpServletRequest.class), any())).thenReturn(testUser);
  }

  @Test
  void testGetMessageOk() throws Exception {

    when(jobsServiceMock.getJob(testUser.getSubjectId(), pipelineId, jobIdOkDone.toString()))
        .thenReturn(jobOkDone);

    mockMvc
        .perform(get(String.format("/api/jobs/v1alpha1/%s/%s", pipelineId, jobIdOkDone)))
        .andExpect(status().isOk())
        .andExpect(content().contentType(MediaType.APPLICATION_JSON));
    //        .andExpect(content().json(jobOkDone.toString())); // this not working
  }

  @Test
  void testGetJobNotFound() throws Exception {
    String badJobId = "bad_job_id";
    when(jobsServiceMock.getJob(testUser.getSubjectId(), pipelineId, badJobId))
        .thenThrow(new JobNotFoundException("some message"));

    mockMvc
        .perform(get(String.format("/api/jobs/v1alpha1/%s/%s", pipelineId, badJobId)))
        .andExpect(status().isNotFound());
  }

  @Test
  void testCreateJobGoodPipeline() throws Exception {
    // This makes the body of the post... which is a lot for very little
    ApiPostJobRequestBody postBody = new ApiPostJobRequestBody().pipelineVersion("TestVersion");
    String postBodyAsJson = MockMvcUtils.convertToJsonString(postBody);

    UUID fakeJobId = UUID.randomUUID();

    // the mocks
    when(jobsServiceMock.createJob(testUser.getSubjectId(), pipelineId, "TestVersion"))
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
    ApiPostJobRequestBody postBody = new ApiPostJobRequestBody().pipelineVersion("TestVersion");
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
    // now that we have the result object, we should further validate the contents of the string
    String contentResults = result.getResponse().getContentAsString();

    // The return value is just the json serialized list of jobs
    ApiGetJobsResponse response =
        new ObjectMapper().readValue(contentResults, ApiGetJobsResponse.class);

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
