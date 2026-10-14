package bio.terra.pipelines.common.utils;

// Note: If there are any changes to this enum, please also update QuotaAllocationSource in
// openapi.yml
public enum QuotaAllocationSourceEnum {
  DEFAULT_FREE,
  FRIENDS_AND_FAMILY_FREE,
  DEVELOPMENT_FREE,
  PRODUCTION_FREE,
  PRODUCTION_CLIENT_TESTING_FREE,
  EXTERNAL_PAID,
  INTERNAL_PAID
}
