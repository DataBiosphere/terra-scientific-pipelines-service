package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.anyString;
import static org.mockito.Mockito.*;

import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.service.exception.InvalidFilterException;
import jakarta.persistence.criteria.*;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.Mock;
import org.mockito.junit.jupiter.MockitoExtension;
import org.springframework.data.jpa.domain.Specification;

@ExtendWith(MockitoExtension.class)
class PipelineRunFilterSpecificationTest {

  @Mock private Root<PipelineRun> root;
  @Mock private CriteriaQuery<?> query;
  @Mock private CriteriaBuilder criteriaBuilder;
  @Mock private Path<Object> path;
  @Mock private Predicate predicate;
  @Mock private Join<Object, Object> pipelineJoin;
  @Mock private Expression<String> stringExpression;

  private static final String TEST_USER_ID = "test-user-123";

  @Test
  void testBuildSpecificationWithUserId_noFilters() {
    when(criteriaBuilder.and(any(Predicate[].class))).thenReturn(predicate);

    Map<String, String> filters = new HashMap<>();
    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, TEST_USER_ID);

    assertNotNull(spec);
    Predicate result = spec.toPredicate(root, query, criteriaBuilder);

    assertNotNull(result);
    verify(criteriaBuilder).equal(root.get("userId"), TEST_USER_ID);
  }

  @Test
  void testBuildSpecificationWithUserId_nullFilters() {
    when(criteriaBuilder.and(any(Predicate[].class))).thenReturn(predicate);

    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(null, TEST_USER_ID);

    assertNotNull(spec);
    Predicate result = spec.toPredicate(root, query, criteriaBuilder);

    assertNotNull(result);
    verify(criteriaBuilder).equal(root.get("userId"), TEST_USER_ID);
  }

  @Test
  void testBuildSpecificationWithUserId_statusFilter_valid() {
    when(criteriaBuilder.and(any(Predicate[].class))).thenReturn(predicate);

    Map<String, String> filters = new HashMap<>();
    filters.put(PipelineRunFilterSpecification.FILTER_STATUS, "RUNNING");

    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, TEST_USER_ID);

    assertNotNull(spec);
    Predicate result = spec.toPredicate(root, query, criteriaBuilder);

    assertNotNull(result);
    verify(criteriaBuilder).equal(root.get("userId"), TEST_USER_ID);
    verify(criteriaBuilder)
        .equal(root.get("status"), CommonPipelineRunStatusEnum.valueOf("RUNNING"));
  }

  @Test
  void testBuildSpecificationWithUserId_statusFilter_caseInsensitive() {
    when(criteriaBuilder.and(any(Predicate[].class))).thenReturn(predicate);

    Map<String, String> filters = new HashMap<>();
    filters.put(PipelineRunFilterSpecification.FILTER_STATUS, "succeeded");

    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, TEST_USER_ID);

    assertNotNull(spec);
    Predicate result = spec.toPredicate(root, query, criteriaBuilder);

    assertNotNull(result);
    verify(criteriaBuilder)
        .equal(root.get("status"), CommonPipelineRunStatusEnum.valueOf("SUCCEEDED"));
  }

  @Test
  void testBuildSpecificationWithUserId_statusFilter_invalid() {
    Map<String, String> filters = new HashMap<>();
    filters.put(PipelineRunFilterSpecification.FILTER_STATUS, "INVALID_STATUS");

    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, TEST_USER_ID);

    InvalidFilterException exception =
        assertThrows(
            InvalidFilterException.class, () -> spec.toPredicate(root, query, criteriaBuilder));

    assertTrue(exception.getMessage().contains("Invalid status. Valid statuses are"));
  }

  @Test
  void testBuildSpecificationWithUserId_jobIdFilter_valid() {
    when(criteriaBuilder.and(any(Predicate[].class))).thenReturn(predicate);

    Map<String, String> filters = new HashMap<>();
    UUID testJobId = UUID.randomUUID();
    filters.put(PipelineRunFilterSpecification.FILTER_JOB_ID, testJobId.toString());

    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, TEST_USER_ID);

    assertNotNull(spec);
    Predicate result = spec.toPredicate(root, query, criteriaBuilder);

    assertNotNull(result);
    verify(criteriaBuilder).equal(root.get("jobId"), testJobId);
  }

  @Test
  void testBuildSpecificationWithUserId_jobIdFilter_invalid() {
    Map<String, String> filters = new HashMap<>();
    filters.put(PipelineRunFilterSpecification.FILTER_JOB_ID, "not-a-uuid");

    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, TEST_USER_ID);

    InvalidFilterException exception =
        assertThrows(
            InvalidFilterException.class, () -> spec.toPredicate(root, query, criteriaBuilder));

    assertTrue(exception.getMessage().contains("Invalid jobId. jobId must be a valid UUID."));
  }

  @Test
  void testBuildSpecificationWithUserId_pipelineNameFilter_valid() {
    when(criteriaBuilder.and(any(Predicate[].class))).thenReturn(predicate);
    when(root.join(anyString())).thenReturn(pipelineJoin);
    when(pipelineJoin.get(anyString())).thenReturn(path);

    Map<String, String> filters = new HashMap<>();
    filters.put(PipelineRunFilterSpecification.FILTER_PIPELINE_NAME, "array_imputation");

    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, TEST_USER_ID);

    assertNotNull(spec);
    Predicate result = spec.toPredicate(root, query, criteriaBuilder);

    assertNotNull(result);
    verify(root).join("pipeline");
    verify(criteriaBuilder).equal(pipelineJoin.get("name"), "array_imputation");
  }

  @Test
  void testBuildSpecificationWithUserId_descriptionFilter() {
    when(criteriaBuilder.and(any(Predicate[].class))).thenReturn(predicate);
    when(criteriaBuilder.lower(any())).thenReturn(stringExpression);

    Map<String, String> filters = new HashMap<>();
    filters.put(PipelineRunFilterSpecification.FILTER_DESCRIPTION, "test description");

    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, TEST_USER_ID);

    assertNotNull(spec);
    Predicate result = spec.toPredicate(root, query, criteriaBuilder);

    assertNotNull(result);
    verify(criteriaBuilder).lower(root.get("description"));
    verify(criteriaBuilder).like(stringExpression, "%test description%");
  }

  @Test
  void testBuildSpecificationWithUserId_multipleFilters() {
    when(criteriaBuilder.and(any(Predicate[].class))).thenReturn(predicate);
    when(criteriaBuilder.lower(any())).thenReturn(stringExpression);

    Map<String, String> filters = new HashMap<>();
    UUID testJobId = UUID.randomUUID();
    filters.put(PipelineRunFilterSpecification.FILTER_STATUS, "RUNNING");
    filters.put(PipelineRunFilterSpecification.FILTER_JOB_ID, testJobId.toString());
    filters.put(PipelineRunFilterSpecification.FILTER_DESCRIPTION, "test");

    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, TEST_USER_ID);

    assertNotNull(spec);
    Predicate result = spec.toPredicate(root, query, criteriaBuilder);

    assertNotNull(result);

    verify(criteriaBuilder).equal(root.get("userId"), TEST_USER_ID);
    verify(criteriaBuilder)
        .equal(root.get("status"), CommonPipelineRunStatusEnum.valueOf("RUNNING"));
    verify(criteriaBuilder).equal(root.get("jobId"), testJobId);
    verify(criteriaBuilder).like(stringExpression, "%test%");
  }

  @Test
  void testBuildSpecificationWithUserId_emptyFilterValue() {
    when(criteriaBuilder.and(any(Predicate[].class))).thenReturn(predicate);

    Map<String, String> filters = new HashMap<>();
    filters.put(PipelineRunFilterSpecification.FILTER_STATUS, "");

    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, TEST_USER_ID);

    assertNotNull(spec);
    Predicate result = spec.toPredicate(root, query, criteriaBuilder);

    assertNotNull(result);

    // Only userId should be filtered, empty value should be ignored
    verify(criteriaBuilder).equal(root.get("userId"), TEST_USER_ID);
  }

  @Test
  void testBuildSpecificationWithUserId_nullFilterValue() {
    when(criteriaBuilder.and(any(Predicate[].class))).thenReturn(predicate);

    Map<String, String> filters = new HashMap<>();
    filters.put(PipelineRunFilterSpecification.FILTER_STATUS, null);

    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, TEST_USER_ID);

    assertNotNull(spec);
    Predicate result = spec.toPredicate(root, query, criteriaBuilder);

    assertNotNull(result);

    // Only userId should be filtered, null value should be ignored
    verify(criteriaBuilder).equal(root.get("userId"), TEST_USER_ID);
  }

  @Test
  void testBuildSpecificationWithUserId_unsupportedFilterKey() {
    Map<String, String> filters = new HashMap<>();
    filters.put("unsupportedKey", "value");

    Specification<PipelineRun> spec =
        PipelineRunFilterSpecification.buildFilterSpecificationWithUserId(filters, TEST_USER_ID);

    InvalidFilterException exception =
        assertThrows(
            InvalidFilterException.class, () -> spec.toPredicate(root, query, criteriaBuilder));

    assertTrue(exception.getMessage().contains("Unsupported filter key"));
    assertTrue(exception.getMessage().contains("unsupportedKey"));
  }
}
