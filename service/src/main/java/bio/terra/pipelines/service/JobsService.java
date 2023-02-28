package bio.terra.pipelines.service;

import bio.terra.pipelines.db.JobsDao;
import bio.terra.pipelines.service.model.Job;
import bio.terra.pipelines.service.model.JobRequest;
import java.sql.Timestamp;
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

  public UUID createJob(String userId, JobRequest jobRequest) {
    UUID jobId = createJobId();
    Timestamp timeSubmitted = getCurrentTimestamp();

    logger.info("Create new {} job with job_id {}", jobRequest.getPipelineId(), jobId);

    // placeholder for actually doing something; for now we're just writing the info to the database
    //        JobBuilder createJob = jobService ...

    String status = "SUBMITTED";
    // Note that this class will grow over time
    // {@link bio/terra/pipelines/service/model/Job.java Job} and {@link bio/terra/pipelines/db/DbJob.java DbJob}
    Job jobFull =
        new Job(
            jobId,
            userId,
            jobRequest.getPipelineId(),
            jobRequest.getPipelineVersion(),
            timeSubmitted,
            null,
            status);

    return jobsDao.createJob(jobFull);
  }

  private UUID createJobId() {
    return UUID.randomUUID();
  }

  private Timestamp getCurrentTimestamp() {
    // TODO add time zone - TSPS-12
    return new Timestamp(System.currentTimeMillis());
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
