package bio.terra.pipelines.common.utils;

import bio.terra.pipelines.db.entities.PipelineRun;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.springframework.data.jpa.domain.Specification;

/** Utility class for building dynamic JPA Specifications for filtering PipelineRun entities */
public class PipelineRunFilterSpecification {

  /**
   * Build a dynamic JPA Specification based on the provided filter map
   *
   * @param filters Map of field names to filter values
   * @param userId The user ID to filter by (always required)
   * @return A Specification for filtering PipelineRun entities
   */
  public static Specification<PipelineRun> buildSpecification(
      Map<String, String> filters, String userId) {
    return (root, query, criteriaBuilder) -> {
      List<jakarta.persistence.criteria.Predicate> predicates = new ArrayList<>();

      // Always filter by userId
      predicates.add(criteriaBuilder.equal(root.get("userId"), userId));

      // Add dynamic filters based on the provided map
      if (filters != null && !filters.isEmpty()) {
        filters.forEach(
            (key, value) -> {
              if (value != null && !value.isEmpty()) {
                switch (key) {
                  case "status":
                    // Convert string to enum
                    try {
                      CommonPipelineRunStatusEnum status =
                          CommonPipelineRunStatusEnum.valueOf(value.toUpperCase());
                      predicates.add(criteriaBuilder.equal(root.get("status"), status));
                    } catch (IllegalArgumentException e) {
                      // Invalid status value, skip this filter
                    }
                    break;
                  case "jobId":
                    // Convert string to UUID
                    try {
                      UUID jobId = UUID.fromString(value);
                      predicates.add(criteriaBuilder.equal(root.get("jobId"), jobId));
                    } catch (IllegalArgumentException e) {
                      // Invalid UUID, skip this filter
                    }
                    break;
                  case "pipelineName":
                    predicates.add(criteriaBuilder.equal(root.get("pipelineName"), value));
                    break;
                  case "description":
                    // Use LIKE for partial matching on description
                    predicates.add(
                        criteriaBuilder.like(
                            criteriaBuilder.lower(root.get("description")),
                            "%" + value.toLowerCase() + "%"));
                    break;
                  default:
                    // Unknown filter field, skip
                    break;
                }
              }
            });
      }

      return criteriaBuilder.and(predicates.toArray(new jakarta.persistence.criteria.Predicate[0]));
    };
  }
}
