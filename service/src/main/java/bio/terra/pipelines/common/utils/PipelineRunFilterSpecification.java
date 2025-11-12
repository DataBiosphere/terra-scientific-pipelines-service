package bio.terra.pipelines.common.utils;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.service.exception.InvalidFilterException;
import jakarta.persistence.criteria.CriteriaBuilder;
import jakarta.persistence.criteria.Join;
import jakarta.persistence.criteria.Predicate;
import jakarta.persistence.criteria.Root;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.springframework.data.jpa.domain.Specification;

/** Utility class for building dynamic JPA Specifications for filtering PipelineRun entities */
public class PipelineRunFilterSpecification {

  public static final String FILTER_STATUS = "status";
  public static final String FILTER_JOB_ID = "jobId";
  public static final String FILTER_PIPELINE_NAME = "pipelineName";
  public static final String FILTER_DESCRIPTION = "description";

  /**
   * Build a dynamic JPA Specification based on the provided map of user-specific filters
   *
   * @param filters Map of field names to filter values
   * @param userId The user ID to filter by (always required)
   * @return A Specification for filtering PipelineRun entities
   */
  public static Specification<PipelineRun> buildFilterSpecificationWithUserId(
      Map<String, String> filters, String userId) {
    return (root, query, criteriaBuilder) -> {
      List<Predicate> predicates = new ArrayList<>();

      // Always filter by userId
      predicates.add(criteriaBuilder.equal(root.get("userId"), userId));

      if (filters != null && !filters.isEmpty()) {
        filters.forEach(
            (key, value) -> {
              if (value != null && !value.isEmpty()) {
                switch (key) {
                  case FILTER_STATUS:
                    predicates.add(validateAndBuildStatusPredicate(value, root, criteriaBuilder));
                    break;
                  case FILTER_JOB_ID:
                    predicates.add(validateAndBuildJobIdPredicate(value, root, criteriaBuilder));
                    break;
                  case FILTER_PIPELINE_NAME:
                    predicates.add(
                        validateAndBuildPipelineNamePredicate(value, root, criteriaBuilder));
                    break;
                  case FILTER_DESCRIPTION:
                    predicates.add(
                        validateAndBuildDescriptionPredicate(value, root, criteriaBuilder));
                    break;
                  default:
                    throw new InvalidFilterException("Unsupported filter key: " + key);
                }
              }
            });
      }

      return criteriaBuilder.and(predicates.toArray(new jakarta.persistence.criteria.Predicate[0]));
    };
  }

  private static Predicate validateAndBuildStatusPredicate(
      String value, Root<PipelineRun> root, CriteriaBuilder criteriaBuilder) {
    try {
      CommonPipelineRunStatusEnum status = CommonPipelineRunStatusEnum.valueOf(value.toUpperCase());
      return criteriaBuilder.equal(root.get(FILTER_STATUS), status);
    } catch (IllegalArgumentException e) {
      throw new InvalidFilterException(
          String.format(
              "Invalid status. Valid statuses are: %s.",
              String.join(
                  ", ",
                  java.util.Arrays.stream(CommonPipelineRunStatusEnum.values())
                      .map(Enum::name)
                      .toArray(String[]::new))));
    }
  }

  private static Predicate validateAndBuildJobIdPredicate(
      String value, Root<PipelineRun> root, CriteriaBuilder criteriaBuilder) {
    try {
      UUID jobId = UUID.fromString(value);
      return criteriaBuilder.equal(root.get(FILTER_JOB_ID), jobId);
    } catch (IllegalArgumentException e) {
      throw new InvalidFilterException("Invalid jobId. jobId must be a valid UUID.");
    }
  }

  private static Predicate validateAndBuildPipelineNamePredicate(
      String value, Root<PipelineRun> root, CriteriaBuilder criteriaBuilder) {
    try {
      PipelinesEnum pipelineName = PipelinesEnum.valueOf(value.toUpperCase());
      // Join to Pipeline table to filter by name
      Join<PipelineRun, Pipeline> pipelineJoin = root.join("pipeline");
      return criteriaBuilder.equal(pipelineJoin.get("name"), pipelineName);
    } catch (IllegalArgumentException e) {
      // this error message intentionally does not list valid pipeline names since some pipelines
      // may be admin-only at times
      throw new InvalidFilterException("Invalid pipeline name filter");
    }
  }

  private static Predicate validateAndBuildDescriptionPredicate(
      String value, Root<PipelineRun> root, CriteriaBuilder criteriaBuilder) {
    return criteriaBuilder.like(
        criteriaBuilder.lower(root.get(FILTER_DESCRIPTION)), "%" + value.toLowerCase() + "%");
  }
}
