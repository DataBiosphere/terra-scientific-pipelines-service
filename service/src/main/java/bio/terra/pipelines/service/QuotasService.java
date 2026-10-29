package bio.terra.pipelines.service;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationSourceEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationStatusEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import bio.terra.pipelines.db.entities.QuotaAllocation;
import bio.terra.pipelines.db.entities.QuotaAllocationDeduction;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.db.repositories.QuotaAllocationDeductionsRepository;
import bio.terra.pipelines.db.repositories.QuotaAllocationsRepository;
import bio.terra.pipelines.db.repositories.UserQuotasRepository;
import bio.terra.pipelines.model.PipelineQuota;
import bio.terra.pipelines.model.QuotaConsumptionRecord;
import bio.terra.pipelines.model.UserQuotaTotals;
import io.sentry.Sentry;
import io.sentry.SentryLevel;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Transactional;

/** This Service handles logic around user quotas and pipeline quotas* */
@Service
public class QuotasService {
  private static final Logger logger = LoggerFactory.getLogger(QuotasService.class);
  private static final String DEFAULT_FREE_UNIQUE_INDEX_NAME =
      "idx_quota_allocations_unique_default_free_per_user_pipeline";
  private final UserQuotasRepository userQuotasRepository;
  private final QuotaAllocationsRepository quotaAllocationsRepository;
  private final QuotaAllocationDeductionsRepository quotaAllocationDeductionsRepository;
  private final PipelineConfigurations pipelineConfigurations;

  @Autowired
  QuotasService(
      UserQuotasRepository userQuotasRepository,
      QuotaAllocationsRepository quotaAllocationsRepository,
      QuotaAllocationDeductionsRepository quotaAllocationDeductionsRepository,
      PipelineConfigurations pipelineConfigurations) {
    this.userQuotasRepository = userQuotasRepository;
    this.quotaAllocationsRepository = quotaAllocationsRepository;
    this.quotaAllocationDeductionsRepository = quotaAllocationDeductionsRepository;
    this.pipelineConfigurations = pipelineConfigurations;
  }

  /** This method gets the PipelineQuota settings for a given pipeline from YAML definition. */
  public PipelineQuota getPipelineQuota(PipelinesEnum pipelineName) {
    PipelineConfigurations.PipelineQuotaConfiguration quotaConfig =
        pipelineConfigurations.getQuotaForPipeline(pipelineName);
    return new PipelineQuota(
        pipelineName,
        quotaConfig.getDefaultQuota(),
        quotaConfig.getMinQuotaConsumed(),
        quotaConfig.getMaxQuotaConsumed(),
        quotaConfig.getQuotaUnits());
  }

  /** This method gets the quota units value for a given pipeline. */
  public QuotaUnitsEnum getQuotaUnitsForPipeline(PipelinesEnum pipelineName) {
    return getPipelineQuota(pipelineName).getQuotaUnits();
  }

  /**
   * This method gets the quota for a given user and pipeline. If the user quota does not exist, it
   * will create a new row in the user quotas table with the default quota for the pipeline.
   */
  public UserQuota getOrCreateQuotaForUserAndPipeline(String userId, PipelinesEnum pipelineName) {
    // try to get the user quota
    Optional<UserQuota> userQuota = getQuotaForUserAndPipeline(userId, pipelineName);

    // if the user quota is not found, grab the default pipeline quota and make a new row in user
    // quotas table
    if (userQuota.isEmpty()) {
      logger.debug(
          "Couldn't find user quota for user {} and pipeline {}. Creating a new row",
          userId,
          pipelineName);
      try {
        /* Note: This try/catch relies on this method NOT being wrapped in an outer @Transactional.
         * Because there is no outer transaction, CrudRepository.save() commits and flushes immediately,
         * allowing us to catch the constraint violation. If you add @Transactional to this service
         * in the future, the flush will be delayed and the try/catch will fail.
         */
        PipelineQuota pipelineQuota = getPipelineQuota(pipelineName);
        UserQuota newUserQuota = new UserQuota();
        newUserQuota.setUserId(userId);
        newUserQuota.setPipelineName(pipelineName);
        newUserQuota.setQuota(pipelineQuota.getDefaultQuota());
        return userQuotasRepository.save(newUserQuota);
      } catch (org.springframework.dao.DataIntegrityViolationException e) {
        // Handle race condition: another thread may have created the quota between our check and
        // insert
        logger.info(
            "Race condition detected while creating user quota for user {} and pipeline {}. Fetching existing record.",
            userId,
            pipelineName);
        return getQuotaForUserAndPipeline(userId, pipelineName)
            .orElseThrow(
                () ->
                    new IllegalStateException(
                        "Quota should exist but was not found for user %s and pipeline %s."
                            .formatted(userId, pipelineName)));
      }
    }
    return userQuota.get();
  }

  /**
   * This method gets the quota for a given user and pipeline. If the user quota does not exist, it
   * will return an empty optional.
   */
  public Optional<UserQuota> getQuotaForUserAndPipeline(String userId, PipelinesEnum pipelineName) {
    return userQuotasRepository.findByUserIdAndPipelineName(userId, pipelineName);
  }

  /**
   * This method updates the quota limit for a given user quota. This should only be called from the
   * Admin Controller
   *
   * @param userQuota - the user quota to update
   * @param newQuotaLimit - the new quota limit
   * @return - the updated user quota
   */
  public UserQuota adminUpdateQuotaLimit(UserQuota userQuota, int newQuotaLimit) {
    if (newQuotaLimit < userQuota.getQuotaConsumed()) {
      throw new InternalServerErrorException(
          String.format(
              "New quota limit: %d, is less than the quota consumed: %d",
              newQuotaLimit, userQuota.getQuotaConsumed()));
    }
    userQuota.setQuota(newQuotaLimit);
    return userQuotasRepository.save(userQuota);
  }

  /**
   * This method validates that the user has enough quota to run the pipeline. If not, it will throw
   * an exception.
   *
   * @param userId - the user id
   * @param pipelineName - the pipeline name
   */
  public void validateUserHasEnoughQuota(String userId, PipelinesEnum pipelineName) {
    UserQuotaTotals userQuotaTotals =
        getOrInitializeQuotaTotalsForUserAndPipeline(userId, pipelineName);
    int minQuotaNeededByPipeline = getPipelineQuota(pipelineName).getMinQuotaConsumed();
    int availableUserQuota = userQuotaTotals.totalAllocated() - userQuotaTotals.totalConsumed();

    if (availableUserQuota < minQuotaNeededByPipeline) {
      throw new BadRequestException(
          ("Insufficient quota to run the pipeline. Quota available: %s, Minimum quota required: %s. "
                  + "Please visit https://services.terra.bio/pipelines/quotas to purchase additional quota.")
              .formatted(availableUserQuota, minQuotaNeededByPipeline));
    }
  }

  /**
   * This method creates a new quota allocation for a given user and pipeline. It will throw an
   * exception if the user already has a DEFAULT_FREE quota allocation for the given pipeline
   *
   * @param userId - the user id
   * @param pipelineName - the pipeline name
   * @param quotaSource - the source this quota allocation should be attributed to
   * @param amount - the amount of quota to allocate
   * @param comments - optional comments to record on the allocation for the audit trail
   * @return - the newly created quota allocation
   */
  public QuotaAllocation allocateQuotaForUserAndPipeline(
      String userId,
      PipelinesEnum pipelineName,
      QuotaAllocationSourceEnum quotaSource,
      int amount,
      String comments) {
    QuotaAllocation newAllocation =
        new QuotaAllocation(
            pipelineName,
            userId,
            quotaSource,
            amount,
            0,
            QuotaAllocationStatusEnum.ACTIVE,
            comments);
    try {
      return quotaAllocationsRepository.save(newAllocation);
    } catch (org.springframework.dao.DataIntegrityViolationException e) {
      if (isDefaultFreeUniqueConstraintViolation(e)) {
        throw new BadRequestException(
            "User already has a DEFAULT_FREE quota allocation for pipeline %s"
                .formatted(pipelineName));
      }
      // If the exception is not DEFAULT_FREE unique index violation, rethrow it
      throw e;
    }
  }

  /**
   * This method returns the total quota allocated and total quota consumed for a given user and
   * pipeline, computed by aggregating across all the user's allocations for that pipeline.
   *
   * @param userId - the user id
   * @param pipelineName - the pipeline name
   * @return - the user's total quota allocated and consumed for the pipeline
   */
  public UserQuotaTotals getQuotaTotalsForUserAndPipeline(
      String userId, PipelinesEnum pipelineName) {
    return quotaAllocationsRepository.sumQuotaTotalsByUserIdAndPipelineName(userId, pipelineName);
  }

  /**
   * This method returns all quota allocations for a given user and pipeline, ordered by created
   * timestamp descending. It will return an empty list if the user has no allocations for the given
   * pipeline.
   *
   * @param userId - the user id
   * @param pipelineName - the pipeline name
   * @return - the list of quota allocations for the user and pipeline
   */
  public List<QuotaAllocation> getAllQuotaAllocationsForUserAndPipeline(
      String userId, PipelinesEnum pipelineName) {
    return quotaAllocationsRepository.findByUserIdAndPipelineNameOrderByCreatedDesc(
        userId, pipelineName);
  }

  /**
   * This method returns the total quota allocated and total quota consumed for a given user and
   * pipeline. If the user has no allocations for the given pipeline AND if the pipeline has a
   * default quota greater than 0, it will create a default allocation for the user/pipeline, and
   * return the totals.
   *
   * @param userId - the user id
   * @param pipelineName - the pipeline name
   * @return - the user's total quota allocated and consumed for the pipeline
   */
  public UserQuotaTotals getOrInitializeQuotaTotalsForUserAndPipeline(
      String userId, PipelinesEnum pipelineName) {
    UserQuotaTotals quotaTotals = getQuotaTotalsForUserAndPipeline(userId, pipelineName);

    // if user doesn’t have any quota AND if pipeline has a default quota greater than 0,
    // create a DEFAULT_FREE allocation
    if (quotaTotals.totalAllocated() == 0) {
      PipelineQuota pipelineQuota = getPipelineQuota(pipelineName);
      if (pipelineQuota.getDefaultQuota() > 0) {
        try {
          allocateQuotaForUserAndPipeline(
              userId,
              pipelineName,
              QuotaAllocationSourceEnum.DEFAULT_FREE,
              pipelineQuota.getDefaultQuota(),
              null);
        } catch (BadRequestException e) {
          // another thread may have created the default allocation between our check and insert
          logger.info(
              "Race condition detected while creating default quota allocation for user {} and pipeline {}. "
                  + "Fetching existing record instead.",
              userId,
              pipelineName);
        }

        // fetch the quota totals again after potentially creating a default allocation
        quotaTotals = getQuotaTotalsForUserAndPipeline(userId, pipelineName);
      }
    }

    return quotaTotals;
  }

  /**
   * This method consumes quota for a given user and pipeline (only from quota_allocations table).
   * It will throw an exception if the user does not have enough quota to consume the requested
   * amount. It will return a list of QuotaConsumptionRecord, one for each allocation that was
   * consumed from.
   *
   * <p>Note: This method is annotated with @Transactional to ensure that all quota consumption
   * operations are atomic and consistent. If any part of the operation fails, the entire
   * transaction will be rolled back, preventing partial consumption of quota.
   *
   * @param userId - the user id
   * @param pipelineName - the pipeline name
   * @param amount - the amount of quota to consume
   * @return - a list of QuotaConsumptionRecord objects
   */
  @Transactional
  public List<QuotaConsumptionRecord> consumeQuotaFromAllocationsForUserAndPipeline(
      String userId, PipelinesEnum pipelineName, int amount) {
    // this should never happen, but if it does, throw an exception
    if (amount <= 0) {
      logger.warn("Received request to consume non-positive amount of quota: {}", amount);
      throw new BadRequestException(
          "Amount to consume must be greater than 0, got %d".formatted(amount));
    }

    // find all active allocations for user/pipeline, ordered by created date (oldest first)
    List<QuotaAllocation> allocations =
        quotaAllocationsRepository.findByUserIdAndPipelineNameAndQuotaStatusOrderByCreatedAsc(
            userId, pipelineName, QuotaAllocationStatusEnum.ACTIVE);

    List<QuotaConsumptionRecord> consumptionRecords = new ArrayList<>();
    int remainingToConsume = amount;

    // do FIFO walk through allocations, consuming quota until the requested amount is consumed
    for (QuotaAllocation allocation : allocations) {
      // if we have consumed the requested amount, break out of the loop
      if (remainingToConsume == 0) {
        break;
      }

      // consume as much as we can from this allocation, up to the remaining amount needed
      // and update the allocation's consumed quota and status as necessary
      int amountConsumedFromAllocation = consumeQuotaFromAllocation(allocation, remainingToConsume);

      // record the consumption in the list of consumption records to return
      consumptionRecords.add(
          new QuotaConsumptionRecord(allocation.getId(), amountConsumedFromAllocation));
      remainingToConsume -= amountConsumedFromAllocation;
    }

    // this should never happen since it should have been verified that user has enough quota for
    // running the job, but if it does, throw an exception
    if (remainingToConsume > 0) {
      throw new BadRequestException(
          "Insufficient quota available to consume for user %s and pipeline %s. Requested: %d, available: %d"
              .formatted(userId, pipelineName, amount, amount - remainingToConsume));
    }

    logger.info(
        "Successfully consumed {} quota for user {} and pipeline {} across {} allocation(s). Allocation(s) consumed: {}",
        amount,
        userId,
        pipelineName,
        consumptionRecords.size(),
        consumptionRecords.stream().map(QuotaConsumptionRecord::quotaAllocationId).toList());

    return consumptionRecords;
  }

  /**
   * This method reverts quota consumption for a list of QuotaConsumptionRecord. It will throw an
   * exception if any of the allocations are not found or if reverting would result in negative
   * consumed quota.
   *
   * <p>It is annotated with @Transactional to ensure that all quota reversion operations are atomic
   * and consistent. If any part of the operation fails, the entire transaction will be rolled back,
   * preventing partial reversion of quota.
   *
   * @param consumptionRecords - a list of QuotaConsumptionRecord objects to revert (must not be
   *     null or empty)
   */
  @Transactional
  public void revertQuotaConsumption(List<QuotaConsumptionRecord> consumptionRecords) {
    // callers are expected to skip the revert when nothing was consumed (which generally
    // should not happen)
    if (consumptionRecords == null || consumptionRecords.isEmpty()) {
      String errorMessage = "No quota consumption records provided to revert.";
      logger.error(errorMessage);
      throw new InternalServerErrorException(errorMessage);
    }

    List<Long> allocationIds =
        consumptionRecords.stream().map(QuotaConsumptionRecord::quotaAllocationId).toList();

    // fetch all allocations together to avoid potential race conditions
    Map<Long, QuotaAllocation> allocationsByIdMap =
        quotaAllocationsRepository.findAllByIdInOrderByCreatedAsc(allocationIds).stream()
            .collect(Collectors.toMap(QuotaAllocation::getId, Function.identity()));

    // revert the quota consumption for each record
    for (QuotaConsumptionRecord consumptionRecord : consumptionRecords) {
      QuotaAllocation allocation = allocationsByIdMap.get(consumptionRecord.quotaAllocationId());

      // this should never happen, but if it does, log an error and throw an exception
      if (allocation == null) {
        throw new IllegalStateException(
            "Quota allocation with ID %d not found while reverting quota consumption."
                .formatted(consumptionRecord.quotaAllocationId()));
      }

      // revert the quota consumption for this allocation
      revertQuotaFromAllocation(allocation, consumptionRecord.amountConsumed());
    }

    // save all the updated allocations to the database in a single batch operation
    quotaAllocationsRepository.saveAll(allocationsByIdMap.values());

    logger.info(
        "Successfully reverted quota consumption for {} allocation(s). Number of allocation(s) reverted: {}",
        consumptionRecords.size(),
        consumptionRecords.stream().map(QuotaConsumptionRecord::quotaAllocationId).toList());
  }

  /**
   * This method records the quota deductions made by a successful pipeline run in the
   * quota_allocation_deductions table, one row per QuotaConsumptionRecord.
   *
   * @param pipelineRunId - the id of the pipeline run that consumed the quota
   * @param consumptionRecords - the records of quota consumed from each allocation
   */
  public void recordQuotaDeductions(
      Long pipelineRunId, List<QuotaConsumptionRecord> consumptionRecords) {
    // a successful run always consumes quota, so missing records indicate a bug upstream
    if (consumptionRecords == null || consumptionRecords.isEmpty()) {
      String errorMessage =
          "No quota consumption records provided to record deductions for pipeline run %d."
              .formatted(pipelineRunId);
      logger.error(errorMessage);
      throw new InternalServerErrorException(errorMessage);
    }

    // this should never happen since we only record deductions for successful runs,
    // but if it does, log a warning and skip recording deductions to avoid duplicate rows
    if (quotaAllocationDeductionsRepository.existsByPipelineRunId(pipelineRunId)) {
      logger.info(
          "Quota deductions already recorded for pipeline run {}. Skipping.", pipelineRunId);
      return;
    }

    // create a list of QuotaAllocationDeduction entities from the consumption records
    List<QuotaAllocationDeduction> deductions =
        consumptionRecords.stream()
            .map(
                consumptionRecord ->
                    new QuotaAllocationDeduction(
                        consumptionRecord.quotaAllocationId(),
                        pipelineRunId,
                        consumptionRecord.amountConsumed(),
                        null))
            .toList();

    // save all the deductions to the database in a single batch operation
    quotaAllocationDeductionsRepository.saveAll(deductions);

    logger.info(
        "Successfully recorded {} quota deduction(s) for pipeline run {}",
        deductions.size(),
        pipelineRunId);
  }

  /**
   * This method consumes quota from a single allocation, up to the amount needed. It updates the
   * allocation's consumed quota and status as necessary, and saves the changes to the database. It
   * returns the amount of quota actually consumed from this allocation.
   *
   * @param allocation - the quota allocation to consume from
   * @param amountNeeded - the amount of quota needed to consume
   * @return - the amount of quota actually consumed from this allocation
   */
  private int consumeQuotaFromAllocation(QuotaAllocation allocation, int amountNeeded) {
    int quotaAvailable = allocation.getQuotaAllocated() - allocation.getQuotaConsumed();

    // this should never happen, but if it does, log an error and throw an exception
    if (quotaAvailable <= 0) {
      String errorMessage =
          "Data integrity issue: allocation %d for user %s pipeline %s is ACTIVE but has 0 or negative quota remaining (%d). Aborting quota consumption."
              .formatted(
                  allocation.getId(),
                  allocation.getUserId(),
                  allocation.getPipelineName(),
                  quotaAvailable);
      logger.error(errorMessage);
      Sentry.captureMessage("Error consuming quota. %s".formatted(errorMessage), SentryLevel.ERROR);
      throw new InternalServerErrorException(errorMessage);
    }

    int amountToConsume = Math.min(quotaAvailable, amountNeeded);
    allocation.setQuotaConsumed(allocation.getQuotaConsumed() + amountToConsume);

    // if there is no remaining quota in that allocation, mark it as EXHAUSTED
    if (allocation.getQuotaConsumed() == allocation.getQuotaAllocated()) {
      allocation.setQuotaStatus(QuotaAllocationStatusEnum.EXHAUSTED);
    }

    quotaAllocationsRepository.save(allocation);

    return amountToConsume;
  }

  /**
   * This method reverts quota consumed from a single allocation by the given amount. It updates the
   * allocation's consumed quota and status as necessary but does not save it; the caller is
   * responsible for persisting the change.
   *
   * @param allocation - the quota allocation to give quota back to
   * @param amountToRevert - the amount of quota to give back to the allocation
   */
  private void revertQuotaFromAllocation(QuotaAllocation allocation, int amountToRevert) {
    int newConsumed = allocation.getQuotaConsumed() - amountToRevert;

    // this should never happen, but if it does, log an error and throw an exception
    if (newConsumed < 0) {
      String errorMessage =
          "Reverting quota consumption would result in negative consumed quota for allocation ID %d."
              .formatted(allocation.getId());
      logger.error(errorMessage);
      throw new IllegalStateException(errorMessage);
    }

    allocation.setQuotaConsumed(newConsumed);

    // if the allocation was previously EXHAUSTED and now has remaining quota, mark it as ACTIVE
    if (allocation.getQuotaStatus() == QuotaAllocationStatusEnum.EXHAUSTED
        && newConsumed < allocation.getQuotaAllocated()) {
      allocation.setQuotaStatus(QuotaAllocationStatusEnum.ACTIVE);
    }
  }

  /**
   * Checks whether the given exception's root cause was thrown by the partial unique index
   * enforcing at most one DEFAULT_FREE allocation per user/pipeline.
   */
  private boolean isDefaultFreeUniqueConstraintViolation(
      org.springframework.dao.DataIntegrityViolationException e) {
    Throwable rootCause = e.getMostSpecificCause();
    return rootCause instanceof org.postgresql.util.PSQLException psqlException
        && psqlException.getServerErrorMessage() != null
        && DEFAULT_FREE_UNIQUE_INDEX_NAME.equals(
            psqlException.getServerErrorMessage().getConstraint());
  }
}
