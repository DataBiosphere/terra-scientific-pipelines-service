package bio.terra.pipelines.service;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.ImputationJob;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.exception.ImputationJobNotFoundException;
import bio.terra.pipelines.db.repositories.ImputationJobsRepository;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoServiceException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.StairwayJobBuilder;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import bio.terra.pipelines.stairway.RunImputationJobFlight;
import com.google.common.annotations.VisibleForTesting;
import java.util.Collections;
import java.util.List;
import java.util.UUID;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
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
  private LeonardoService leonardoService;
  private SamService samService;
  private WdsService wdsService;
  private final StairwayJobService stairwayJobService;
  private ImputationConfiguration imputationConfiguration;

  @Autowired
  ImputationService(
      ImputationJobsRepository imputationJobsRepository,
      PipelineInputsRepository pipelineInputsRepository,
      LeonardoService leonardoService,
      SamService samService,
      WdsService wdsService,
      StairwayJobService stairwayJobService,
      ImputationConfiguration imputationConfiguration) {
    this.imputationJobsRepository = imputationJobsRepository;
    this.pipelineInputsRepository = pipelineInputsRepository;
    this.leonardoService = leonardoService;
    this.samService = samService;
    this.wdsService = wdsService;
    this.stairwayJobService = stairwayJobService;
    this.imputationConfiguration = imputationConfiguration;
  }

  /**
   * Creates a new Imputation pipeline service job, using a Stairway flight, based on a user's
   * request. Returns jobId of new job (which is the same as the flightId) if flight submission is
   * successful.
   *
   * @param userId
   * @param pipelineVersion
   * @return String jobId
   *     <p>Note that the information in the requested job will grow over time, along with the
   *     following related classes:
   * @see ImputationJob
   */
  @Transactional
  public UUID createImputationJob(String userId, String pipelineVersion, Object pipelineInputs) {
    logger.info("Create new imputation version {} job for user {}", pipelineVersion, userId);

    StairwayJobBuilder stairwayJobBuilder =
        stairwayJobService
            .newJob()
            .jobId(createJobId())
            .flightClass(RunImputationJobFlight.class)
            .pipelineId(PipelinesEnum.IMPUTATION)
            .pipelineVersion(pipelineVersion)
            .userId(userId)
            .pipelineInputs(pipelineInputs);

    return stairwayJobBuilder.submit();
  }

  @Transactional
  public ImputationJob getImputationJob(UUID jobId, String userId) {
    return imputationJobsRepository
        .findJobByJobIdAndUserId(jobId, userId)
        .orElseThrow(
            () ->
                new ImputationJobNotFoundException(
                    String.format("ImputationJob %s for user %s not found.", jobId, userId)));
  }

  @VisibleForTesting
  protected List<ImputationJob> getImputationJobs(String userId) {
    return imputationJobsRepository.findAllByUserId(userId);
  }

  // TSPS-136 will require that the user provide the job UUID, and should remove this method
  protected UUID createJobId() {
    return UUID.randomUUID();
  }

  @Transactional
  public UUID writeJobToDb(
      UUID jobUuid, String userId, String pipelineVersion, Object pipelineInputs) {

    // write job to imputation database
    ImputationJob job = new ImputationJob();
    job.setJobId(jobUuid);
    job.setUserId(userId);
    job.setPipelineVersion(pipelineVersion);

    ImputationJob createdJob = writeJobToDbThrowsDuplicateException(job);

    // save related pipeline inputs
    PipelineInput pipelineInput = new PipelineInput();
    pipelineInput.setJobId(createdJob.getId());
    pipelineInput.setInputs(pipelineInputs.toString());
    pipelineInputsRepository.save(pipelineInput);

    return createdJob.getJobId();
  }

  public List<ListAppResponse> queryForWorkspaceApps() {
    String workspaceId = imputationConfiguration.workspaceId();
    try {
      List<ListAppResponse> getAppsResponse =
          leonardoService.getApps(workspaceId, samService.getTspsServiceAccountToken(), false);

      logger.info(
          "GetAppsResponse for workspace id {}: {}",
          imputationConfiguration.workspaceId(),
          getAppsResponse);

      String wdsUri =
          leonardoService.getWdsUrlFromApps(
              workspaceId, samService.getTspsServiceAccountToken(), false);
      logger.info("Wds uri for workspace id {}: {}", imputationConfiguration.workspaceId(), wdsUri);
      logger.info(
          "Wds health: {}",
          wdsService.checkHealth(wdsUri, samService.getTspsServiceAccountToken()));

      logger.info(
          "Wds schema: {}",
          wdsService.querySchema(
              wdsUri,
              samService.getTspsServiceAccountToken(),
              imputationConfiguration.workspaceId()));

      return getAppsResponse;
    } catch (LeonardoServiceException e) {
      logger.error("Get Apps called for workspace id {} failed", workspaceId);
      return Collections.emptyList();
    } catch (WdsServiceException e) {
      logger.error("Calls to Wds for workspace id {} failed", workspaceId);
      return Collections.emptyList();
    }
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
