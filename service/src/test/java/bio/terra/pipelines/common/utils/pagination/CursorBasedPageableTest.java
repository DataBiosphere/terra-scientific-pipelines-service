package bio.terra.pipelines.common.utils.pagination;

import static bio.terra.pipelines.common.utils.pagination.CursorBasedPageable.ENCODED_PADDING_AROUND_VALUE;
import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.testutils.BaseTest;
import java.time.LocalDateTime;
import java.util.Base64;
import java.util.Objects;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

class CursorBasedPageableTest extends BaseTest {

  private static final String ENCODED_TEN_VALUE =
      "IyMjMTAjIyMgLSAyMDI0LTA3LTIzVDE0OjQzOjEyLjQ5NDk4OA==";

  @Test
  void decodeValue() {
    String originalValue = "10";
    CursorBasedPageable cursorBasedPageable = new CursorBasedPageable(100, null, null);
    // encode value
    String encodedString = CursorBasedPageable.getEncodedCursor(originalValue, true);
    // decode value
    String decodedValue = CursorBasedPageable.getDecodedCursor(encodedString);
    assertEquals(originalValue, decodedValue);
  }

  @Test
  void decodeNullValueOrEmpty() {
    CursorBasedPageable cursorBasedPageable = new CursorBasedPageable(100, null, null);
    // decode value
    assertThrows(IllegalArgumentException.class, () -> CursorBasedPageable.getDecodedCursor(null));
    assertThrows(IllegalArgumentException.class, () -> CursorBasedPageable.getDecodedCursor(""));
    assertThrows(
        IllegalArgumentException.class, () -> CursorBasedPageable.getDecodedCursor("    "));
  }

  @Test
  void encodeValue() {
    CursorBasedPageable cursorBasedPageable = new CursorBasedPageable(10, null, null);

    // because this uses a LocalDateTime.now() as part of its encoding strategy, its hard to test
    // the full string matches so we are just testing the first 10 characters as a "good enough"
    String structuredValue =
        ENCODED_PADDING_AROUND_VALUE
            + 10
            + ENCODED_PADDING_AROUND_VALUE
            + " - "
            + LocalDateTime.now();
    String encodedString = Base64.getEncoder().encodeToString(structuredValue.getBytes());
    assertEquals(
        encodedString.substring(0, 10),
        Objects.requireNonNull(CursorBasedPageable.getEncodedCursor("10", true)).substring(0, 10));
  }

  @Test
  void encodeValueNoPreviousOrNextPages() {
    CursorBasedPageable cursorBasedPageable = new CursorBasedPageable(10, null, null);
    assertFalse(cursorBasedPageable.hasCursors());
    assertNull(CursorBasedPageable.getEncodedCursor("10", false));
  }

  @Test
  void hasNextCursor() {
    CursorBasedPageable cursorBasedPageableNoNext = new CursorBasedPageable(10, null, "blah");
    Assertions.assertFalse(cursorBasedPageableNoNext.hasNextPageCursor());
    assertTrue(cursorBasedPageableNoNext.hasCursors());
    CursorBasedPageable cursorBasedPageableNext = new CursorBasedPageable(10, "balhblah", null);
    assertTrue(cursorBasedPageableNext.hasNextPageCursor());
    assertEquals("balhblah", cursorBasedPageableNext.getNextPageCursor());
  }

  @Test
  void hasPrevCursor() {
    CursorBasedPageable cursorBasedPageableNoPrev = new CursorBasedPageable(10, "blah", null);
    Assertions.assertFalse(cursorBasedPageableNoPrev.hasPrevPageCursor());
    assertTrue(cursorBasedPageableNoPrev.hasCursors());
    CursorBasedPageable cursorBasedPageablePrev = new CursorBasedPageable(10, null, "balhblah");
    assertTrue(cursorBasedPageablePrev.hasPrevPageCursor());
    assertEquals("balhblah", cursorBasedPageablePrev.getPrevPageCursor());
  }

  @Test
  void getSearchValue() {
    CursorBasedPageable cursorBasedPageableNoCursor = new CursorBasedPageable(10, null, null);
    assertEquals(
        CursorBasedPageable.POSTGRES_INT4_MAX_VALUE, cursorBasedPageableNoCursor.getSearchValue());
    CursorBasedPageable cursorBasedPageable = new CursorBasedPageable(10, ENCODED_TEN_VALUE, null);
    assertEquals("10", cursorBasedPageable.getSearchValue());
    cursorBasedPageable = new CursorBasedPageable(10, null, ENCODED_TEN_VALUE);
    assertEquals("10", cursorBasedPageable.getSearchValue());
  }

  @Test
  void cursorBasedPageableHashCode() {
    CursorBasedPageable cursorBasedPageableNextCursor =
        new CursorBasedPageable(10, "hmmm_ok", null);
    // 17 and 31 are hardcoded in this hashCode method of this class
    assertEquals(
        new HashCodeBuilder(17, 31).append(10).append("hmmm_ok").append((String) null).toHashCode(),
        cursorBasedPageableNextCursor.hashCode());

    CursorBasedPageable cursorBasedPageablePrevCursor =
        new CursorBasedPageable(10, null, "hmmm_ok");
    assertEquals(
        new HashCodeBuilder(17, 31).append(10).append((String) null).append("hmmm_ok").toHashCode(),
        cursorBasedPageablePrevCursor.hashCode());
  }

  @Test
  void cursorBasedPageableEquals() {
    CursorBasedPageable cursorBasedPageableOne = new CursorBasedPageable(10, "hmmm_ok", null);
    CursorBasedPageable cursorBasedPageableSameAsOne = new CursorBasedPageable(10, "hmmm_ok", null);
    CursorBasedPageable cursorBasedPageableDifferentFromOne =
        new CursorBasedPageable(15, null, "hmmm_ok");

    assertEquals(cursorBasedPageableOne, cursorBasedPageableSameAsOne);
    assertNotEquals(cursorBasedPageableOne, cursorBasedPageableDifferentFromOne);
    assertNotEquals(cursorBasedPageableSameAsOne, cursorBasedPageableDifferentFromOne);
  }
}
