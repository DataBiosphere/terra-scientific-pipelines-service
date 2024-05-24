package bio.terra.pipelines.service;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Job;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.JobsRepository;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlight;
import bio.terra.pipelines.stairway.imputation.RunImputationJobFlightMapKeys;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.util.Map;
import java.util.UUID;
import java.util.stream.Collectors;
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
  private final JobsRepository jobsRepository;
  private final PipelineInputsRepository pipelineInputsRepository;
  private final JobService jobService;

  private final ObjectMapper objectMapper = new ObjectMapper();

  @Autowired
  ImputationService(
      JobsRepository jobsRepository,
      PipelineInputsRepository pipelineInputsRepository,
      JobService jobService) {
    this.jobsRepository = jobsRepository;
    this.pipelineInputsRepository = pipelineInputsRepository;
    this.jobService = jobService;
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
   * @see Job
   */
  public UUID createImputationJob(
      UUID jobId,
      String userId,
      String description,
      Pipeline imputationPipeline,
      Map<String, Object> userProvidedPipelineInputs,
      String resultPath) {

    PipelinesEnum imputationPipelineName =
        PipelinesEnum.valueOf(imputationPipeline.getName().toUpperCase());
    logger.info("Create new {} job for user {}", imputationPipelineName, userId);

    JobBuilder jobBuilder =
        jobService
            .newJob()
            .jobId(jobId)
            .flightClass(RunImputationJobFlight.class)
            .addParameter(JobMapKeys.PIPELINE_NAME.getKeyName(), imputationPipelineName)
            .addParameter(JobMapKeys.USER_ID.getKeyName(), userId)
            .addParameter(JobMapKeys.DESCRIPTION.getKeyName(), description)
            .addParameter(RunImputationJobFlightMapKeys.PIPELINE_ID, imputationPipeline.getId())
            .addParameter(
                RunImputationJobFlightMapKeys.PIPELINE_INPUT_DEFINITIONS,
                // we need to do this to get the object to deserialize correctly
                imputationPipeline.getPipelineInputDefinitions().stream()
                    .collect(Collectors.toList()))
            .addParameter(
                RunImputationJobFlightMapKeys.USER_PROVIDED_PIPELINE_INPUTS,
                userProvidedPipelineInputs)
            .addParameter(
                RunImputationJobFlightMapKeys.CONTROL_WORKSPACE_ID,
                imputationPipeline.getWorkspaceId().toString())
            .addParameter(
                RunImputationJobFlightMapKeys.WDL_METHOD_NAME,
                imputationPipeline.getWdlMethodName())
            .addParameter(JobMapKeys.RESULT_PATH.getKeyName(), resultPath);

    return jobBuilder.submit();
  }

  @Transactional
  public UUID writeJobToDb(
      UUID jobUuid, String userId, Long pipelineId, Map<String, Object> pipelineInputs) {

    // write job to database
    Job job = new Job();
    job.setJobId(jobUuid);
    job.setUserId(userId);
    job.setPipelineId(pipelineId);

    Job createdJob = writeJobToDbThrowsDuplicateException(job);

    String pipelineInputsAsString;
    try {
      // do this to write the pipeline inputs without writing the class name
      pipelineInputsAsString = objectMapper.writeValueAsString(pipelineInputs);
    } catch (JsonProcessingException e) {
      throw new InternalServerErrorException(
          "THIS SHOULD NEVER HAPPEN! Error converting pipeline inputs to string", e);
    }

    // save related pipeline inputs
    PipelineInput pipelineInput = new PipelineInput();
    pipelineInput.setJobId(createdJob.getId());
    pipelineInput.setInputs(pipelineInputsAsString);
    pipelineInputsRepository.save(pipelineInput);

    return createdJob.getJobId();
  }

  protected Job writeJobToDbThrowsDuplicateException(Job job) throws DuplicateObjectException {
    try {
      jobsRepository.save(job);
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
