package bio.terra.pipelines.controller;

import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.when;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.common.iam.BearerToken;
import bio.terra.common.iam.BearerTokenFactory;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.controller.JobsApiController;
import bio.terra.pipelines.config.SamConfiguration;
import bio.terra.pipelines.db.exception.JobNotFoundException;
import bio.terra.pipelines.iam.SamService;
import bio.terra.pipelines.service.JobsService;
import bio.terra.pipelines.service.model.Job;
import java.sql.Timestamp;
import java.util.Optional;
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

@ContextConfiguration(classes = JobsApiController.class)
@WebMvcTest
class JobsControllerTest {
  @MockBean JobsService serviceMock;
  @MockBean SamUserFactory samUserFactoryMock;
  @MockBean BearerTokenFactory bearerTokenFactory;
  @MockBean SamConfiguration samConfiguration;
  @MockBean SamService samService;

  @Autowired private MockMvc mockMvc;
  private SamUser testUser =
      new SamUser(
          "test@email",
          UUID.randomUUID().toString(),
          new BearerToken(UUID.randomUUID().toString()));

  private String pipelineId = "imputation";
  private String badJobId = "bad_job_id";
  private Timestamp timestamp = new Timestamp(System.currentTimeMillis());
  private UUID jobIdOkSubmitted = UUID.randomUUID();
  private Job jobOkSubmitted =
      new Job(
          jobIdOkSubmitted,
          testUser.getSubjectId(),
          pipelineId,
          "v0",
          timestamp,
          Optional.empty(),
          "Submitted");

  private UUID jobIdOkDone = UUID.randomUUID();
  private Job jobOkDone =
      new Job(
          jobIdOkDone,
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

    when(serviceMock.getJob(testUser.getSubjectId(), pipelineId, jobIdOkDone.toString()))
        .thenReturn(jobOkDone);

    mockMvc
        .perform(get(String.format("/api/jobs/v1alpha1/%s/%s", pipelineId, jobIdOkDone)))
        .andExpect(status().isOk())
        .andExpect(content().contentType(MediaType.APPLICATION_JSON));
    //        .andExpect(content().json(jobOkDone.toString())); // this not working
  }

  @Test
  void testGetJobNotFound() throws Exception {
    when(serviceMock.getJob(testUser.getSubjectId(), pipelineId, badJobId))
        .thenThrow(new JobNotFoundException("some message"));

    mockMvc
        .perform(get(String.format("/api/jobs/v1alpha1/%s/%s", pipelineId, badJobId)))
        .andExpect(status().isNotFound());
  }
}
