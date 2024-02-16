package bio.terra.pipelines.service;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.ImputationJob;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.ImputationJobsRepository;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlight;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import java.util.UUID;
import org.hibernate.exception.ConstraintViolationException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DataIntegrityViolationException;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Transactional;

/** Service to encapsulate logic used to run an imputation pipeline */
@Service
public class ImputationService {
  private static final Logger logger = LoggerFactory.getLogger(ImputationService.class);
  private final ImputationJobsRepository imputationJobsRepository;
  private final PipelineInputsRepository pipelineInputsRepository;
  private final JobService jobService;
  private ImputationConfiguration imputationConfiguration;

  @Autowired
  ImputationService(
      ImputationJobsRepository imputationJobsRepository,
      PipelineInputsRepository pipelineInputsRepository,
      JobService jobService,
      ImputationConfiguration imputationConfiguration) {
    this.imputationJobsRepository = imputationJobsRepository;
    this.pipelineInputsRepository = pipelineInputsRepository;
    this.jobService = jobService;
    this.imputationConfiguration = imputationConfiguration;
  }

  /**
   * Creates a new Imputation pipeline service job, using a Stairway flight, based on a user's
   * request. Returns jobId of new job (which is the same as the flightId) if flight submission is
   * successful.
   *
   * @param userId
   * @param imputationPipeline - a pipeline that handles imputation
   * @return String jobId
   *     <p>Note that the information in the requested job will grow over time, along with the
   *     following related classes:
   * @see ImputationJob
   */
  public UUID createImputationJob(
      UUID jobId,
      String userId,
      String description,
      Pipeline imputationPipeline,
      Object pipelineInputs,
      String resultPath) {

    PipelinesEnum pipelineName = PipelinesEnum.valueOf(imputationPipeline.getName().toUpperCase());
    logger.info("Create new {} job for user {}", pipelineName, userId);

    JobBuilder jobBuilder =
        jobService
            .newJob()
            .jobId(jobId)
            .flightClass(RunImputationJobFlight.class)
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), pipelineName)
            .addParameter(JobMapKeys.USER_ID.getKeyName(), userId)
            .addParameter(JobMapKeys.DESCRIPTION.getKeyName(), description)
            .addParameter(RunImputationJobFlightMapKeys.PIPELINE_ID, imputationPipeline.getId())
            .addParameter(RunImputationJobFlightMapKeys.PIPELINE_INPUTS, pipelineInputs)
            .addParameter(
                RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID,
                imputationConfiguration.getWorkspaceId())
            .addParameter(
                RunImputationJobFlightMapKeys.WDL_METHOD_NAME,
                imputationPipeline.getWdlMethodName())
            .addParameter(JobMapKeys.RESULT_PATH.getKeyName(), resultPath);

    return jobBuilder.submit();
  }

  @Transactional
  public UUID writeJobToDb(UUID jobUuid, String userId, Long pipelineId, Object pipelineInputs) {

    // write job to imputation database
    ImputationJob job = new ImputationJob();
    job.setJobId(jobUuid);
    job.setUserId(userId);
    job.setPipelineId(pipelineId);

    ImputationJob createdJob = writeJobToDbThrowsDuplicateException(job);

    // save related pipeline inputs
    PipelineInput pipelineInput = new PipelineInput();
    pipelineInput.setJobId(createdJob.getId());
    pipelineInput.setInputs(pipelineInputs.toString());
    pipelineInputsRepository.save(pipelineInput);

    return createdJob.getJobId();
  }

  protected ImputationJob writeJobToDbThrowsDuplicateException(ImputationJob job)
      throws DuplicateObjectException {
    try {
      imputationJobsRepository.save(job);
      logger.info("job saved for jobId: {}", job.getJobId());
    } catch (DataIntegrityViolationException e) {
      if (e.getCause() instanceof ConstraintViolationException c
          && c.getConstraintName().contains("jobId_unique")) {
        throw new DuplicateObjectException(
            String.format("Duplicate jobId %s found", job.getJobId()));
      }
      throw e;
    }

    return job;
  }
}
