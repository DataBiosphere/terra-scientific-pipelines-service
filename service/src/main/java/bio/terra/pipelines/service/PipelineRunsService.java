package bio.terra.pipelines.service;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.db.entities.PipelineInput;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.repositories.PipelineInputsRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.util.Map;
import java.util.UUID;
import org.hibernate.exception.ConstraintViolationException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DataIntegrityViolationException;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Transactional;

/** Service to encapsulate logic used to manage pipeline runs */
@Service
public class PipelineRunsService {
  private static final Logger logger = LoggerFactory.getLogger(PipelineRunsService.class);
  private final PipelineRunsRepository pipelineRunsRepository;
  private final PipelineInputsRepository pipelineInputsRepository;

  private final ObjectMapper objectMapper = new ObjectMapper();

  @Autowired
  public PipelineRunsService(
      PipelineRunsRepository pipelineRunsRepository,
      PipelineInputsRepository pipelineInputsRepository) {
    this.pipelineRunsRepository = pipelineRunsRepository;
    this.pipelineInputsRepository = pipelineInputsRepository;
  }

  @Transactional
  public UUID writePipelineRunToDb(
      UUID jobUuid, String userId, Long pipelineId, Map<String, Object> pipelineInputs) {

    // write pipelineRun to database
    PipelineRun pipelineRun = new PipelineRun();
    pipelineRun.setJobId(jobUuid);
    pipelineRun.setUserId(userId);
    pipelineRun.setPipelineId(pipelineId);

    PipelineRun createdPipelineRun = writePipelineRunToDbThrowsDuplicateException(pipelineRun);

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
    pipelineInput.setJobId(createdPipelineRun.getId());
    pipelineInput.setInputs(pipelineInputsAsString);
    pipelineInputsRepository.save(pipelineInput);

    return createdPipelineRun.getJobId();
  }

  protected PipelineRun writePipelineRunToDbThrowsDuplicateException(PipelineRun pipelineRun)
      throws DuplicateObjectException {
    try {
      pipelineRunsRepository.save(pipelineRun);
      logger.info("pipelineRun saved for jobId: {}", pipelineRun.getJobId());
    } catch (DataIntegrityViolationException e) {
      if (e.getCause() instanceof ConstraintViolationException c
          && c.getConstraintName().contains("jobId_unique")) {
        throw new DuplicateObjectException(
            String.format("Duplicate jobId %s found", pipelineRun.getJobId()));
      }
      throw e;
    }

    return pipelineRun;
  }
}
