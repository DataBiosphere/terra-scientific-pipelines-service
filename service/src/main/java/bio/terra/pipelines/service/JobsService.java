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

  @Autowired
  public JobsService(JobsDao jobsDao) {
    this.jobsDao = jobsDao;
  }

  /**
   * Creates a new pipeline service job based on a user's request. Returns jobId of new job if write
   * to db is successful, otherwise returns null.
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
    Instant timeSubmitted = getCurrentTimestamp();

    logger.info("Create new {} version {} job for user {}", pipelineId, pipelineVersion, userId);

    // placeholder for actually doing something; for now we're just writing the info to the database
    //        JobBuilder createJob = jobService ...

    String status = "SUBMITTED";

    return writeJobToDb(userId, pipelineId, pipelineVersion, timeSubmitted, status, 1);
  }

  protected UUID createJobId() {
    return UUID.randomUUID();
  }

  protected UUID writeJobToDb(
      String userId,
      String pipelineId,
      String pipelineVersion,
      Instant timeSubmitted,
      String status,
      int attempt) {

    UUID jobUuid = createJobId();

    Job jobToStore =
        new Job(jobUuid, userId, pipelineId, pipelineVersion, timeSubmitted, null, status);

    if (attempt <= 3) {
      UUID createdJobId = jobsDao.createJob(jobToStore);
      if (createdJobId == null) {
        int nextAttempt = attempt + 1;
        return writeJobToDb(
            userId, pipelineId, pipelineVersion, timeSubmitted, status, nextAttempt);
      } else {
        return createdJobId;
      }
    } else {
      // 3 attempts to write a job to the database failed
      return null;
    }
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
