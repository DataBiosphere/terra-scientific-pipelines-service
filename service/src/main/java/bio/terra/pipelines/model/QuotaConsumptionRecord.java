package bio.terra.pipelines.model;

/**
 * Represents a record of quota consumption from a quota allocation.
 *
 * @param quotaAllocationId - the ID of the quota allocation
 * @param amountConsumed - the amount of quota consumed in this record
 */
public record QuotaConsumptionRecord(Long quotaAllocationId, int amountConsumed) {}
