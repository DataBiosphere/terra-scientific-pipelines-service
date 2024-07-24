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
 * This class is used to help build a predicate to both filter by greater than or less than a field
 * name and value as well as ordering the results of the query, currently this is hardcoded to
 * descending to get the most recent records first. This is done based on the CursorBasePageable
 * passed in.
 *
 * @param <T> - DB entity you are filtering on
 */
@RequiredArgsConstructor
public class PageSpecification<T> implements Specification<T> {

  private final transient String mainFieldName;
  private final transient CursorBasedPageable cursorBasedPageable;

  @Override
  public Predicate toPredicate(
      @NotNull Root<T> root, CriteriaQuery<?> query, @NotNull CriteriaBuilder criteriaBuilder) {
    var predicate = applyPaginationFilter(root, criteriaBuilder);
    query.orderBy(criteriaBuilder.desc(root.get(mainFieldName)));

    return predicate;
  }

  private Predicate applyPaginationFilter(Root<T> root, CriteriaBuilder criteriaBuilder) {
    var searchValue = cursorBasedPageable.getSearchValue();

    return cursorBasedPageable.hasPrevPageCursor()
        ? criteriaBuilder.greaterThan(root.get(mainFieldName), searchValue)
        : criteriaBuilder.lessThan(root.get(mainFieldName), searchValue);
  }

  public int hashCode() {
    return new HashCodeBuilder(17, 31)
        // two randomly chosen prime numbers
        // if deriving: appendSuper(super.hashCode()).
        .append(mainFieldName)
        .append(cursorBasedPageable)
        .toHashCode();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof PageSpecification<?> otherObject)) return false;
    if (obj == this) return true;

    return new EqualsBuilder()
        .append(mainFieldName, otherObject.mainFieldName)
        .append(cursorBasedPageable, otherObject.cursorBasedPageable)
        .isEquals();
  }
}
