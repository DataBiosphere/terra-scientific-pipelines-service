package bio.terra.pipelines.common.utils;

import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.service.exception.InvalidFilterException;
import jakarta.persistence.criteria.Join;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.springframework.data.jpa.domain.Specification;

/** Utility class for building dynamic JPA Specifications for filtering PipelineRun entities */
public class PipelineRunFilterSpecification {

  /**
   * Build a dynamic JPA Specification based on the provided map of user-specific filters
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

      if (filters != null && !filters.isEmpty()) {
        filters.forEach(
            (key, value) -> {
              if (value != null && !value.isEmpty()) {
                switch (key) {
                  case "status":
                    try {
                      CommonPipelineRunStatusEnum status =
                          CommonPipelineRunStatusEnum.valueOf(value.toUpperCase());
                      predicates.add(criteriaBuilder.equal(root.get("status"), status));
                    } catch (IllegalArgumentException e) {
                      throw new InvalidFilterException(
                          "Invalid status value. Valid statuses are: "
                              + String.join(
                                  ", ",
                                  java.util.Arrays.stream(CommonPipelineRunStatusEnum.values())
                                      .map(Enum::name)
                                      .toArray(String[]::new))
                              + ". "
                              + value);
                    }
                    break;
                  case "jobId":
                    try {
                      UUID jobId = UUID.fromString(value);
                      predicates.add(criteriaBuilder.equal(root.get("jobId"), jobId));
                    } catch (IllegalArgumentException e) {
                      throw new InvalidFilterException(
                          "Invalid jobId format. jobId must be a UUID." + value);
                    }
                    break;
                  case "pipelineName":
                    try {
                      PipelinesEnum pipelineName = PipelinesEnum.valueOf(value.toUpperCase());
                      // Join to Pipeline table to filter by name
                      Join<PipelineRun, Pipeline> pipelineJoin = root.join("pipeline");
                      predicates.add(criteriaBuilder.equal(pipelineJoin.get("name"), pipelineName));
                    } catch (InvalidFilterException e) {
                      throw new IllegalArgumentException("Invalid pipeline name: " + value);
                    }
                    break;
                  case "description":
                    predicates.add(
                        criteriaBuilder.like(
                            criteriaBuilder.lower(root.get("description")),
                            "%" + value.toLowerCase() + "%"));
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
}
