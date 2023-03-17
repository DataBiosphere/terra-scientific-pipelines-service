package bio.terra.pipelines.service;

import bio.terra.pipelines.db.JobsDao;
import bio.terra.pipelines.service.model.Job;
import java.time.Instant;
import java.util.List;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

/** The Jobs Service manages job requests to run the service's Scientific Pipelines. */
@Component
public class JobsService {

  private static final Logger logger = LoggerFactory.getLogger(JobsService.class);

  private final JobsDao jobsDao;
  private final PipelinesService pipelinesService;

  @Autowired
  public JobsService(JobsDao jobsDao, PipelinesService pipelinesService) {
    this.jobsDao = jobsDao;
    this.pipelinesService = pipelinesService;
  }

  /**
   * Creates a new pipeline service job based on a user's request.
   *
   * @param userId
   * @param pipelineId
   * @param pipelineVersion
   * @return UUID jobId
   *     <p>Note that the information in the requested job will grow over time, along with the
   *     following related classes:
   * @see bio.terra.pipelines.db.JobsDao
   * @see bio.terra.pipelines.service.model.Job
   */
  public UUID createJob(String userId, String pipelineId, String pipelineVersion) {
    // validate that the requested pipelineId exists
    pipelinesService.validatePipeline(pipelineId);

    UUID jobId = createJobId();
    Instant timeSubmitted = getCurrentTimestamp();

    logger.info("Create new {} version {} job with job_id {}", pipelineId, pipelineVersion, jobId);

    // placeholder for actually doing something; for now we're just writing the info to the database
    //        JobBuilder createJob = jobService ...

    String status = "SUBMITTED";
    Job jobFull = new Job(jobId, userId, pipelineId, pipelineVersion, timeSubmitted, null, status);

    return jobsDao.createJob(jobFull);
  }

  protected UUID createJobId() {
    return UUID.randomUUID();
  }

  private Instant getCurrentTimestamp() {
    // Instant creates a timestamp in UTC
    return Instant.now();
  }

  public List<Job> getJobs(String userId, String pipelineId) {
    logger.info("Get all jobs in {} pipeline for user {}}", pipelineId, userId);
    return jobsDao.getJobs(userId, pipelineId);
  }

  public Job getJob(String userId, String pipelineId, String jobId) {
    logger.info("Get job {} in {} pipeline for user {}}", jobId, pipelineId, userId);
    return jobsDao.getJob(userId, pipelineId, jobId);
  }
}
