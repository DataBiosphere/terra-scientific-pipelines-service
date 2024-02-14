package bio.terra.pipelines.service;

import bio.terra.cbas.model.*;
import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.ImputationJob;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.ImputationJobsRepository;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.cbas.CbasServiceException;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoServiceException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import bio.terra.pipelines.stairway.RunImputationJobFlight;
import bio.terra.pipelines.stairway.RunImputationJobFlightMapKeys;
import java.util.Collections;
import java.util.List;
import java.util.UUID;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
import org.databiosphere.workspacedata.model.RecordAttributes;
import org.databiosphere.workspacedata.model.RecordRequest;
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
  private CbasService cbasService;
  private final JobService jobService;
  private ImputationConfiguration imputationConfiguration;

  @Autowired
  ImputationService(
      ImputationJobsRepository imputationJobsRepository,
      PipelineInputsRepository pipelineInputsRepository,
      LeonardoService leonardoService,
      SamService samService,
      WdsService wdsService,
      CbasService cbasService,
      JobService jobService,
      ImputationConfiguration imputationConfiguration) {
    this.imputationJobsRepository = imputationJobsRepository;
    this.pipelineInputsRepository = pipelineInputsRepository;
    this.leonardoService = leonardoService;
    this.samService = samService;
    this.wdsService = wdsService;
    this.cbasService = cbasService;
    this.jobService = jobService;
    this.imputationConfiguration = imputationConfiguration;
  }

  /**
   * Creates a new Imputation pipeline service job, using a Stairway flight, based on a user's
   * request. Returns jobId of new job (which is the same as the flightId) if flight submission is
   * successful.
   *
   * @param userId
   * @param imputationPipeline - a pipeline that handlines imputation
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

  public List<ListAppResponse> queryForWorkspaceApps(Pipeline pipeline, UUID flightUUID) {
    String workspaceId = imputationConfiguration.workspaceId();
    try {
      // leonardo related calls
      List<ListAppResponse> getAppsResponse =
          leonardoService.getApps(workspaceId, samService.getTspsServiceAccountToken(), false);

      logger.info(
          "GetAppsResponse for workspace id {}: {}",
          imputationConfiguration.workspaceId(),
          getAppsResponse);

      //      // example to create cromwell_runner_app
      //      CreateAppRequest createAppRequest =
      //          new CreateAppRequest()
      //              .appType(AppType.CROMWELL_RUNNER_APP)
      //              .accessScope(AppAccessScope.USER_PRIVATE)
      //              .labels(Map.of("saturnAutoCreated", "true"));
      //      logger.info(
      //          "creating runner app for workspace {}: {}",
      //          workspaceId,
      //          leonardoService.createAppV2(
      //              workspaceId,
      //              samService.getTspsServiceAccountToken(),
      //              "tsps-cr-" + UUID.randomUUID(),
      //              createAppRequest));

      // wds related calls
      String wdsUri =
          leonardoService.getWdsUrlFromApps(
              workspaceId, samService.getTspsServiceAccountToken(), false);

      logger.info("Wds uri for workspace id {}: {}", imputationConfiguration.workspaceId(), wdsUri);

      logger.info(
          "Wds health: {}",
          wdsService.checkHealth(wdsUri, samService.getTspsServiceAccountToken()));

      String id = UUID.randomUUID().toString();
      String wdsTableName = "hello_world";
      String primaryKey = "workflow_id";
      RecordAttributes recordAttributes = new RecordAttributes();
      recordAttributes.put("scatter", 2);
      RecordRequest createRecordRequest = new RecordRequest().attributes(recordAttributes);
      logger.info(
          "add record to {} table: {}",
          wdsTableName,
          wdsService.createOrReplaceRecord(
              wdsUri,
              samService.getTspsServiceAccountToken(),
              createRecordRequest,
              imputationConfiguration.workspaceId(),
              wdsTableName,
              id,
              primaryKey));

      logger.info(
          "Wds schema: {}",
          wdsService.querySchema(
              wdsUri,
              samService.getTspsServiceAccountToken(),
              imputationConfiguration.workspaceId()));

      // cbas related calls
      String cbasUri =
          leonardoService.getCbasUrlFromApps(
              workspaceId, samService.getTspsServiceAccountToken(), false);

      logger.info(
          "cbas uri for workspace id {}: {}", imputationConfiguration.workspaceId(), cbasUri);

      logger.info(
          "Cbas health: {}",
          cbasService.checkHealth(cbasUri, samService.getTspsServiceAccountToken()));

      //      // example of creating method using cbas service
      //      PostMethodRequest postMethodRequest =
      //          new PostMethodRequest()
      //              .methodName(pipeline.getWdlMethodName())
      //              .methodDescription("method description")
      //              .methodSource(PostMethodRequest.MethodSourceEnum.GITHUB)
      //              .methodUrl(pipeline.getWdlUrl())
      //              .methodVersion("1.0");
      //      logger.info(
      //          "this is creating a new method in cbas: {}",
      //          cbasService.createMethod(
      //              cbasUri, samService.getTspsServiceAccountToken(), postMethodRequest));

      MethodListResponse methodListResponse =
          cbasService.getAllMethods(cbasUri, samService.getTspsServiceAccountToken());
      logger.info("list of methods available: {}", methodListResponse);

      // grab methodVersionId needed to submit a submission
      UUID methodVersionId = null;
      for (MethodDetails methodDetails : methodListResponse.getMethods()) {
        if (methodDetails.getName().equals(pipeline.getWdlMethodName())) {
          methodVersionId = methodDetails.getMethodVersions().get(0).getMethodVersionId();
          logger.info("this is the method version id: {}", methodVersionId);
          break;
        }
      }

      // launch a cbas submission
      RunSetRequest runSetRequest =
          new RunSetRequest()
              .runSetDescription(String.format("flight id: %s", flightUUID))
              .runSetName("run set name")
              .methodVersionId(methodVersionId)
              .addWorkflowInputDefinitionsItem(
                  new WorkflowInputDefinition()
                      .inputName("HelloWorld.scatter_num")
                      .inputType(
                          new ParameterTypeDefinitionPrimitive()
                              .primitiveType(PrimitiveParameterValueType.INT)
                              .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                      .source(
                          new ParameterDefinitionRecordLookup()
                              .recordAttribute("scatter")
                              .type(ParameterDefinition.TypeEnum.RECORD_LOOKUP)))
              .addWorkflowOutputDefinitionsItem(
                  new WorkflowOutputDefinition()
                      .outputName("HelloWorld.output_file")
                      .outputType(
                          new ParameterTypeDefinitionArray()
                              .nonEmpty(true)
                              .arrayType(
                                  new ParameterTypeDefinitionPrimitive()
                                      .primitiveType(PrimitiveParameterValueType.FILE)
                                      .type(ParameterTypeDefinition.TypeEnum.PRIMITIVE))
                              .type(ParameterTypeDefinition.TypeEnum.ARRAY))
                      .destination(
                          new OutputDestinationRecordUpdate()
                              .recordAttribute("output_file")
                              .type(OutputDestination.TypeEnum.RECORD_UPDATE)))
              .wdsRecords(new WdsRecordSet().recordType(wdsTableName).addRecordIdsItem(id));
      logger.info(
          "run set created: {}",
          cbasService.createRunSet(
              cbasUri, samService.getTspsServiceAccountToken(), runSetRequest));

      return getAppsResponse;
    } catch (LeonardoServiceException e) {
      logger.error("Get Apps called for workspace id {} failed", workspaceId);
      return Collections.emptyList();
    } catch (WdsServiceException e) {
      logger.error("Calls to Wds for workspace id {} failed", workspaceId);
      return Collections.emptyList();
    } catch (CbasServiceException e) {
      logger.error("Calls to Cbas for workspace id {} failed", workspaceId);
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
