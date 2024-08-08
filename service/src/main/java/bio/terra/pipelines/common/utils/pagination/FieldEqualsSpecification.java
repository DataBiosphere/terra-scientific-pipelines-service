package bio.terra.pipelines.common.utils.pagination;

import jakarta.persistence.criteria.CriteriaBuilder;
import jakarta.persistence.criteria.CriteriaQuery;
import jakarta.persistence.criteria.Predicate;
import jakarta.persistence.criteria.Root;
import lombok.RequiredArgsConstructor;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.jetbrains.annotations.NotNull;
import org.springframework.data.jpa.domain.Specification;

/**
 * This class is used to help build a predicate to filter based on a match to a value for a column.
 * This is currently used to filter paginated results based on the user id field but can be used
 * generically for any column / value combination.
 *
 * @param <T> - DB entity you are filtering on
 */
@RequiredArgsConstructor
public class FieldEqualsSpecification<T> implements Specification<T> {

  private final transient String fieldName;
  private final transient String fieldValue;

  @Override
  public Predicate toPredicate(
      @NotNull Root<T> root, CriteriaQuery<?> query, @NotNull CriteriaBuilder criteriaBuilder) {
    return criteriaBuilder.equal(root.get(fieldName), fieldValue);
  }

  public int hashCode() {
    return new HashCodeBuilder(17, 31)
        // two randomly chosen prime numbers
        // if deriving: appendSuper(super.hashCode()).
        .append(fieldName)
        .append(fieldValue)
        .toHashCode();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof FieldEqualsSpecification<?> otherObject)) return false;
    if (obj == this) return true;

    return new EqualsBuilder()
        .append(fieldName, otherObject.fieldName)
        .append(fieldValue, otherObject.fieldValue)
        .isEquals();
  }
}
