package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.ApiException;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.config.SamConfiguration;
import bio.terra.pipelines.db.entities.Job;
import bio.terra.pipelines.db.exception.PipelineNotFoundException;
import bio.terra.pipelines.generated.api.JobsApi;
import bio.terra.pipelines.generated.model.*;
import bio.terra.pipelines.service.JobsService;
import bio.terra.pipelines.service.PipelinesService;
import io.swagger.annotations.Api;
import java.util.List;
import java.util.UUID;
import javax.servlet.http.HttpServletRequest;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestBody;

/** Jobs controller */
@Controller
@Api(tags = {"jobs"})
public class JobsApiController implements JobsApi {
  private final SamConfiguration samConfiguration;
  private final SamUserFactory samUserFactory;
  private final HttpServletRequest request;
  private final JobsService jobsService;
  private final PipelinesService pipelinesService;

  @Autowired
  public JobsApiController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      JobsService jobsService,
      PipelinesService pipelinesService) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.jobsService = jobsService;
    this.pipelinesService = pipelinesService;
  }

  private static final Logger logger = LoggerFactory.getLogger(JobsApiController.class);

  private SamUser getAuthenticatedInfo() {
    return samUserFactory.from(request, samConfiguration.basePath());
  }

  // -- Jobs --

  @Override
  public ResponseEntity<ApiPostJobResponse> createJob(
      @PathVariable("pipelineId") String pipelineId, @RequestBody ApiPostJobRequestBody body) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    String pipelineVersion = body.getPipelineVersion();
    Object pipelineInputs = body.getPipelineInputs();

    if (!pipelinesService.pipelineExists(pipelineId)) {
      throw new PipelineNotFoundException(
          String.format("Requested pipeline %s not found.", pipelineId));
    }

    logger.info(
        "Creating {} pipeline job (version {}) for {} subject {} with inputs {}",
        pipelineId,
        pipelineVersion,
        userRequest.getEmail(),
        userId,
        pipelineInputs);

    // TODO assuming we will write outputs back to source workspace, we will need to check user
    // permissions for write access to the workspace - explore interceptors

    UUID createdJobUuid = jobsService.createJob(userId, pipelineId, pipelineVersion);
    if (createdJobUuid == null) {
      logger.error("New {} pipeline job creation failed.", pipelineId);
      throw new ApiException("An internal error occurred.");
    }

    ApiPostJobResponse createdJobResponse = new ApiPostJobResponse();
    createdJobResponse.setJobId(createdJobUuid.toString());
    logger.info("Created {} job {}", pipelineId, createdJobUuid);

    return new ResponseEntity<>(createdJobResponse, HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiGetJobResponse> getJob(
      @PathVariable("pipelineId") String pipelineId, @PathVariable("jobId") UUID jobId) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    Job job = jobsService.getJob(userId, pipelineId, jobId);
    ApiGetJobResponse result = jobToApi(job);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiGetJobsResponse> getJobs(@PathVariable("pipelineId") String pipelineId) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    List<Job> jobList = jobsService.getJobs(userId, pipelineId);
    ApiGetJobsResponse result = jobsToApi(jobList);

    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  static ApiGetJobResponse jobToApi(Job job) {
    ApiGetJobResponse apiGetJobResponse =
        new ApiGetJobResponse()
            .jobId(job.getJobId().toString())
            .userId(job.getUserId())
            .pipelineId(job.getPipelineId())
            .pipelineVersion(job.getPipelineVersion())
            .timeSubmitted(job.getTimeSubmitted().toString())
            .status(job.getStatus());
    if (job.getTimeCompleted() != null) {
      apiGetJobResponse.setTimeCompleted(job.getTimeCompleted().toString());
    }
    return apiGetJobResponse;
  }

  static ApiGetJobsResponse jobsToApi(List<Job> jobList) {
    ApiGetJobsResponse apiResult = new ApiGetJobsResponse();

    for (Job job : jobList) {
      var apiJob = jobToApi(job);

      apiResult.add(apiJob);
    }

    return apiResult;
  }
}
