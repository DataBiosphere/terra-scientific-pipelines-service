package bio.terra.pipelines.service;

import bio.terra.pipelines.db.entities.Job;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.exception.JobNotFoundException;
import bio.terra.pipelines.db.repositories.JobsRepository;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.dependencies.stairway.StairwayJobBuilder;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.stairway.CreateJobFlight;
import java.time.Instant;
import java.util.List;
import java.util.UUID;
import org.hibernate.exception.ConstraintViolationException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DataIntegrityViolationException;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Transactional;

/** The Jobs Service manages job requests to run the service's Scientific Pipelines. */
@Service
public class JobsService {

  private static final Logger logger = LoggerFactory.getLogger(JobsService.class);
  private final JobsRepository jobsRepository;
  private final PipelineInputsRepository pipelineInputsRepository;
  private final StairwayJobService stairwayJobService;

  @Autowired
  public JobsService(
      JobsRepository jobsRepository,
      PipelineInputsRepository pipelineInputsRepository,
      StairwayJobService stairwayJobService) {
    this.jobsRepository = jobsRepository;
    this.pipelineInputsRepository = pipelineInputsRepository;
    this.stairwayJobService = stairwayJobService;
  }

  /**
   * Creates a new pipeline service job, using a Stairway flight, based on a user's request. Returns
   * jobId of new job (which is the same as the flightId) if flight submission is successful.
   *
   * @param userId
   * @param pipelineId
   * @param pipelineVersion
   * @return String jobId
   *     <p>Note that the information in the requested job will grow over time, along with the
   *     following related classes:
   * @see Job
   */
  @Transactional
  public String createJob(
      String userId, String pipelineId, String pipelineVersion, Object pipelineInputs) {
    logger.info("Create new {} version {} job for user {}", pipelineId, pipelineVersion, userId);

    StairwayJobBuilder stairwayJobBuilder =
        stairwayJobService
            .newJob()
            .jobId(createJobId())
            .flightClass(CreateJobFlight.class)
            .pipelineId(pipelineId)
            .pipelineVersion(pipelineVersion)
            .submittingUserId(userId)
            .pipelineInputs(pipelineInputs);

    return stairwayJobBuilder.submit();
  }

  protected String createJobId() {
    return UUID.randomUUID().toString();
  }

  public UUID writeJobToDb(
      UUID jobUuid,
      String userId,
      String pipelineId,
      String pipelineVersion,
      Instant timeSubmitted,
      String status,
      Object pipelineInputs) {

    Job job = new Job();
    job.setJobId(jobUuid);
    job.setUserId(userId);
    job.setPipelineId(pipelineId);
    job.setPipelineVersion(pipelineVersion);
    job.setTimeSubmitted(timeSubmitted);
    job.setTimeCompleted(null);
    job.setStatus(status);

    Job createdJob = writeJobToDbThrowsDuplicateException(job);

    // once job is created save related pipeline inputs
    PipelineInput pipelineInput = new PipelineInput();
    pipelineInput.setJobId(createdJob.getId());
    pipelineInput.setInputs(pipelineInputs.toString());
    pipelineInputsRepository.save(pipelineInput);

    return createdJob.getJobId();
  }

  protected Job writeJobToDbThrowsDuplicateException(Job job) throws DuplicateObjectException {
    try {
      jobsRepository.save(job);
      logger.info("job saved for jobId: {}", job.getJobId());
    } catch (DataIntegrityViolationException e) {
      if (e.getCause() instanceof ConstraintViolationException) {
        throw new DuplicateObjectException(
            String.format("Duplicate jobId %s found", job.getJobId()));
      }
      throw e;
    }

    return job;
  }

  public Instant getCurrentTimestamp() {
    // Instant creates a timestamp in UTC
    return Instant.now();
  }

  public List<Job> getJobs(String userId, String pipelineId) {
    logger.info("Get all jobs in {} pipeline for user {}}", pipelineId, userId);
    return jobsRepository.findAllByPipelineIdAndUserId(pipelineId, userId);
  }

  public Job getJob(String userId, String pipelineId, UUID jobId) {
    logger.info("Get job {} in {} pipeline for user {}}", jobId, pipelineId, userId);
    return jobsRepository
        .findJobByPipelineIdAndUserIdAndJobId(pipelineId, userId, jobId)
        .orElseThrow(() -> new JobNotFoundException(String.format("Job %s not found.", jobId)));
  }
}
