package bio.terra.pipelines.common.utils.pagination;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

import bio.terra.pipelines.testutils.BaseTest;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.junit.jupiter.api.Test;

class FieldEqualsSpecificationTest extends BaseTest {

  @Test
  void fieldEqualsSpecificationHashCode() {
    FieldEqualsSpecification<?> fieldEqualsSpecification =
        new FieldEqualsSpecification<>("userId", "12");
    assertEquals(
        new HashCodeBuilder(17, 31).append("userId").append("12").toHashCode(),
        fieldEqualsSpecification.hashCode());
  }

  @Test
  void fieldEqualsSpecificationEquals() {
    FieldEqualsSpecification<?> fieldEqualsSpecificationOne =
        new FieldEqualsSpecification<>("userId", "12");
    FieldEqualsSpecification<?> fieldEqualsSpecificationSameAsOne =
        new FieldEqualsSpecification<>("userId", "12");
    FieldEqualsSpecification<?> fieldEqualsSpecificationDifferentFromOne =
        new FieldEqualsSpecification<>("userId", "43");

    assertEquals(fieldEqualsSpecificationOne, fieldEqualsSpecificationSameAsOne);
    assertNotEquals(fieldEqualsSpecificationOne, fieldEqualsSpecificationDifferentFromOne);
    assertNotEquals(fieldEqualsSpecificationSameAsOne, fieldEqualsSpecificationDifferentFromOne);
  }
}
