package bio.terra.pipelines.service;

import bio.terra.common.exception.NotFoundException;
import bio.terra.pipelines.common.utils.DataDeliveryStatusEnum;
import bio.terra.pipelines.db.entities.DataDelivery;
import bio.terra.pipelines.db.repositories.DataDeliveryRepository;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

/** Service to manage data delivery records for pipeline runs */
@Service
public class DataDeliveryService {
  private static final Logger logger = LoggerFactory.getLogger(DataDeliveryService.class);
  private final DataDeliveryRepository dataDeliveryRepository;

  @Autowired
  public DataDeliveryService(DataDeliveryRepository dataDeliveryRepository) {
    this.dataDeliveryRepository = dataDeliveryRepository;
  }

  /**
   * Create and save a new data delivery record
   *
   * @param pipelineRunId - the pipeline run ID
   * @param jobId - the job ID
   * @param status - the delivery status
   * @param gcsDestinationPath - the GCS destination path
   * @return the saved DataDelivery entity
   */
  public DataDelivery createDataDelivery(
      Long pipelineRunId, UUID jobId, DataDeliveryStatusEnum status, String gcsDestinationPath) {
    logger.info(
        "Created data delivery record for pipeline run ID {} and job ID {}", pipelineRunId, jobId);

    DataDelivery dataDelivery = new DataDelivery(pipelineRunId, jobId, status, gcsDestinationPath);
    return dataDeliveryRepository.save(dataDelivery);
  }

  /**
   * Get the latest data delivery record for a specific pipeline run, ordered by creation time
   *
   * @param pipelineRunId - the pipeline run ID
   * @return the most recent DataDelivery entity, or null if none exist
   */
  public DataDelivery getLatestDataDeliveryByPipelineRunId(Long pipelineRunId) {
    return dataDeliveryRepository
        .findFirstByPipelineRunIdOrderByCreatedDesc(pipelineRunId)
        .orElse(null);
  }

  /**
   * Update the status of an existing data delivery record
   *
   * @param pipelineRunId - the job ID of the record to update
   * @param newStatus - the new status value
   * @return the updated DataDelivery entity
   * @throws NotFoundException if no record exists for the job ID
   */
  public DataDelivery updateDataDeliveryStatus(
      Long pipelineRunId, DataDeliveryStatusEnum newStatus) {
    DataDelivery dataDelivery = getLatestDataDeliveryByPipelineRunId(pipelineRunId);

    logger.info(
        "Updating data delivery record with ID {} for pipeline run ID {} from status {} to status {}",
        dataDelivery != null ? dataDelivery.getId() : "null",
        pipelineRunId,
        dataDelivery != null ? dataDelivery.getStatus() : "null",
        newStatus);

    if (dataDelivery == null) {
      String errorMessage =
          String.format("No data delivery record found for pipeline run ID %s", pipelineRunId);
      logger.error(errorMessage);
      throw new NotFoundException(errorMessage);
    }

    dataDelivery.setStatus(newStatus);
    return dataDeliveryRepository.save(dataDelivery);
  }
}
