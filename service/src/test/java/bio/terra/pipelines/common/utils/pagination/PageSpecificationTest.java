package bio.terra.pipelines.common.utils.pagination;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

import bio.terra.pipelines.testutils.BaseTest;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.junit.jupiter.api.Test;

class PageSpecificationTest extends BaseTest {

  @Test
  void pageSpecificationHashCode() {
    CursorBasedPageable cursorBasedPageable = new CursorBasedPageable(14, "doesnt_matter", null);
    PageSpecification<?> pageSpecification = new PageSpecification<>("id", cursorBasedPageable);
    assertEquals(
        new HashCodeBuilder(17, 31).append("id").append(cursorBasedPageable).toHashCode(),
        pageSpecification.hashCode());
  }

  @Test
  void pageSpecificationEquals() {
    CursorBasedPageable cursorBasedPageable = new CursorBasedPageable(14, "doesnt_matter", null);
    PageSpecification<?> pageSpecificationOne = new PageSpecification<>("id", cursorBasedPageable);
    PageSpecification<?> pageSpecificationSameAsOne =
        new PageSpecification<>("id", cursorBasedPageable);
    PageSpecification<?> pageSpecificationDifferentFromOne =
        new PageSpecification<>("not_id", cursorBasedPageable);

    assertEquals(pageSpecificationOne, pageSpecificationSameAsOne);
    assertNotEquals(pageSpecificationOne, pageSpecificationDifferentFromOne);
    assertNotEquals(pageSpecificationSameAsOne, pageSpecificationDifferentFromOne);
  }
}
