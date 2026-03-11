package bio.terra.pipelines.service;

import bio.terra.common.exception.NotFoundException;
import bio.terra.pipelines.db.entities.DataDelivery;
import bio.terra.pipelines.db.repositories.DataDeliveryRepository;
import java.util.List;
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
      Long pipelineRunId, UUID jobId, String status, String gcsDestinationPath) {
    logger.info(
        "Creating data delivery record for pipelineRunId: {}, jobId: {}", pipelineRunId, jobId);

    DataDelivery dataDelivery = new DataDelivery(pipelineRunId, jobId, status, gcsDestinationPath);
    DataDelivery savedDelivery = dataDeliveryRepository.save(dataDelivery);

    logger.info("Data delivery record created with id: {}", savedDelivery.getId());
    return savedDelivery;
  }

  /**
   * Get a data delivery record by job ID
   *
   * @param jobId - the job ID to search for
   * @return the DataDelivery entity if found
   * @throws NotFoundException if no record exists for the job ID
   */
  public DataDelivery getDataDeliveryByJobId(UUID jobId) {
    logger.info("Retrieving data delivery record for jobId: {}", jobId);
    return dataDeliveryRepository
        .findByJobId(jobId)
        .orElseThrow(
            () ->
                new NotFoundException(
                    String.format("No data delivery record found for jobId: %s", jobId)));
  }

  /**
   * Get all data delivery records for a specific pipeline run
   *
   * @param pipelineRunId - the pipeline run ID
   * @return list of DataDelivery entities for the pipeline run
   */
  public List<DataDelivery> getDataDeliveriesByPipelineRunId(Long pipelineRunId) {
    logger.info("Retrieving data delivery records for pipelineRunId: {}", pipelineRunId);
    return dataDeliveryRepository.findAllByPipelineRunId(pipelineRunId);
  }

  /**
   * Get the latest data delivery record for a specific pipeline run, ordered by creation time
   *
   * @param pipelineRunId - the pipeline run ID
   * @return the most recent DataDelivery entity, or null if none exist
   */
  public DataDelivery getLatestDataDeliveryByPipelineRunId(Long pipelineRunId) {
    logger.info("Retrieving latest data delivery record for pipelineRunId: {}", pipelineRunId);
    return dataDeliveryRepository
        .findFirstByPipelineRunIdOrderByCreatedDesc(pipelineRunId)
        .orElse(null);
  }

  /**
   * Get all data delivery records with a specific status
   *
   * @param status - the status to filter by
   * @return list of DataDelivery entities with the specified status
   */
  public List<DataDelivery> getDataDeliveriesByStatus(String status) {
    logger.info("Retrieving data delivery records with status: {}", status);
    return dataDeliveryRepository.findAllByStatus(status);
  }

  /**
   * Update the status of an existing data delivery record
   *
   * @param pipelineRunId - the job ID of the record to update
   * @param newStatus - the new status value
   * @return the updated DataDelivery entity
   * @throws NotFoundException if no record exists for the job ID
   */
  public DataDelivery updateDataDeliveryStatus(Long pipelineRunId, String newStatus) {
    logger.info(
        "Updating data delivery status for pipelineRunId: {} to {}", pipelineRunId, newStatus);

    DataDelivery dataDelivery =
        getDataDeliveriesByPipelineRunId(pipelineRunId).stream().findFirst().orElse(null);
    dataDelivery.setStatus(newStatus);
    DataDelivery updatedDelivery = dataDeliveryRepository.save(dataDelivery);

    logger.info("Data delivery status updated for pipelineRunId: {}", pipelineRunId);
    return updatedDelivery;
  }

  /**
   * Update the GCS destination path of an existing data delivery record
   *
   * @param jobId - the job ID of the record to update
   * @param newGcsDestinationPath - the new GCS destination path
   * @return the updated DataDelivery entity
   * @throws NotFoundException if no record exists for the job ID
   */
  public DataDelivery updateDataDeliveryGcsPath(UUID jobId, String newGcsDestinationPath) {
    logger.info("Updating GCS destination path for jobId: {} to {}", jobId, newGcsDestinationPath);

    DataDelivery dataDelivery = getDataDeliveryByJobId(jobId);
    dataDelivery.setGcsDestinationPath(newGcsDestinationPath);
    DataDelivery updatedDelivery = dataDeliveryRepository.save(dataDelivery);

    logger.info("GCS destination path updated for jobId: {}", jobId);
    return updatedDelivery;
  }

  /**
   * Update both status and GCS destination path of an existing data delivery record
   *
   * @param jobId - the job ID of the record to update
   * @param newStatus - the new status value
   * @param newGcsDestinationPath - the new GCS destination path
   * @return the updated DataDelivery entity
   * @throws NotFoundException if no record exists for the job ID
   */
  public DataDelivery updateDataDelivery(
      UUID jobId, String newStatus, String newGcsDestinationPath) {
    logger.info(
        "Updating data delivery for jobId: {} - status: {}, path: {}",
        jobId,
        newStatus,
        newGcsDestinationPath);

    DataDelivery dataDelivery = getDataDeliveryByJobId(jobId);
    dataDelivery.setStatus(newStatus);
    dataDelivery.setGcsDestinationPath(newGcsDestinationPath);
    DataDelivery updatedDelivery = dataDeliveryRepository.save(dataDelivery);

    logger.info("Data delivery updated for jobId: {}", jobId);
    return updatedDelivery;
  }

  /**
   * Check if a data delivery record exists for a job ID
   *
   * @param jobId - the job ID to check
   * @return true if a record exists, false otherwise
   */
  public boolean existsByJobId(UUID jobId) {
    return dataDeliveryRepository.existsByJobId(jobId);
  }
}
