package bio.terra.pipelines.model;

/**
 * Aggregate of a user's total quota allocated and consumed across all of their allocations for a
 * given pipeline.
 *
 * @param totalAllocated - sum of quota_allocated across all user's allocations for a given pipeline
 * @param totalConsumed - sum of quota_consumed across all user's allocations for a given pipeline
 */
public record UserQuotaTotals(int totalAllocated, int totalConsumed) {}
