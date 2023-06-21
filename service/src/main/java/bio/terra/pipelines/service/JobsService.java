package bio.terra.pipelines.service;

import bio.terra.pipelines.db.entities.DbJob;
import bio.terra.pipelines.db.exception.JobNotFoundException;
import bio.terra.pipelines.db.repositories.JobsRepository;
import java.time.Instant;
import java.util.List;
import java.util.Optional;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

/** The Jobs Service manages job requests to run the service's Scientific Pipelines. */
@Service
public class JobsService {

  private static final Logger logger = LoggerFactory.getLogger(JobsService.class);
  private final JobsRepository jobsRepository;

  @Autowired
  public JobsService(JobsRepository jobsRepository) {
    this.jobsRepository = jobsRepository;
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
   * @see bio.terra.pipelines.db.entities.DbJob
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

    DbJob dbJob = new DbJob();
    dbJob.setJobId(jobUuid.toString());
    dbJob.setUserId(userId);
    dbJob.setPipelineId(pipelineId);
    dbJob.setPipelineVersion(pipelineVersion);
    dbJob.setTimeSubmitted(timeSubmitted);
    dbJob.setTimeCompleted(null);
    dbJob.setStatus(status);

    if (attempt <= 3) {
      UUID createdJobId = writeJobToDbRetryDuplicateException(dbJob);
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

  protected UUID writeJobToDbRetryDuplicateException(DbJob dbJob) {
    Optional<DbJob> jobExists = jobsRepository.findJobByJobId(dbJob.getJobId());
    if (jobExists.isPresent()) {
      logger.warn("Duplicate jobId {} found, retrying", dbJob.getJobId());
      return null;
    }
    jobsRepository.save(dbJob);
    logger.info("job saved for jobId: {}", dbJob.getJobId());

    return UUID.fromString(dbJob.getJobId());
  }

  private Instant getCurrentTimestamp() {
    // Instant creates a timestamp in UTC
    return Instant.now();
  }

  public List<DbJob> getJobs(String userId, String pipelineId) {
    logger.info("Get all jobs in {} pipeline for user {}}", pipelineId, userId);
    return jobsRepository.findAllByPipelineIdAndUserId(pipelineId, userId);
  }

  public DbJob getJob(String userId, String pipelineId, String jobId) {
    logger.info("Get job {} in {} pipeline for user {}}", jobId, pipelineId, userId);
    return jobsRepository
        .findJobByPipelineIdAndUserIdAndJobId(pipelineId, userId, jobId)
        .orElseThrow(() -> new JobNotFoundException(String.format("Job %s not found.", jobId)));
  }
}
