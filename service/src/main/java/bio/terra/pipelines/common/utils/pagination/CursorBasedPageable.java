package bio.terra.pipelines.common.utils.pagination;

import static java.util.Objects.requireNonNull;
import static org.apache.commons.lang3.StringUtils.substringBetween;

import java.time.LocalDateTime;
import java.util.Base64;
import lombok.AllArgsConstructor;
import lombok.Data;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;

/**
 * This class is used to help with determining page token encoding/decoding for our paginated
 * endpoints. There is logic to decode already created page tokens and encode the field you are
 * paginating over into page tokens. It is assumed that the field you are paginating over is a
 * postgres id field of type int4
 */
@Data
@AllArgsConstructor
public class CursorBasedPageable {
  private int size;
  private final String nextPageCursor;
  private final String prevPageCursor;

  public static final String POSTGRES_INT4_MAX_VALUE = "2147483647";

  public boolean hasNextPageCursor() {
    return nextPageCursor != null && !nextPageCursor.isEmpty();
  }

  public boolean hasPrevPageCursor() {
    return prevPageCursor != null && !prevPageCursor.isEmpty();
  }

  public boolean hasCursors() {
    return hasPrevPageCursor() || hasNextPageCursor();
  }

  public String getDecodedCursor(String cursorValue) {
    if (cursorValue == null || cursorValue.isEmpty()) {
      throw new IllegalArgumentException("Cursor value is not valid!");
    }
    var decodedBytes = Base64.getDecoder().decode(cursorValue);
    var decodedValue = new String(decodedBytes);

    return substringBetween(decodedValue, "###");
  }

  public String getEncodedCursor(String field, boolean hasPrevOrNextElements) {
    requireNonNull(field);

    if (!hasPrevOrNextElements) return null;

    var structuredValue = "###" + field + "### - " + LocalDateTime.now();
    return Base64.getEncoder().encodeToString(structuredValue.getBytes());
  }

  public String getSearchValue() {
    // postgres does not have a positive infinity symbol for an int4 field type so we hard code the
    // largest possible
    // value that can occur in that column type here.  very hacky
    if (!hasCursors()) return POSTGRES_INT4_MAX_VALUE;

    return hasPrevPageCursor()
        ? getDecodedCursor(prevPageCursor)
        : getDecodedCursor(nextPageCursor);
  }

  @Override
  public int hashCode() {
    return new HashCodeBuilder(17, 31)
        // two randomly chosen prime numbers
        // if deriving: appendSuper(super.hashCode()).
        .append(size)
        .append(nextPageCursor)
        .append(prevPageCursor)
        .toHashCode();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof CursorBasedPageable otherObject)) return false;
    if (obj == this) return true;

    return new EqualsBuilder()
        .append(size, otherObject.size)
        .append(nextPageCursor, otherObject.nextPageCursor)
        .append(prevPageCursor, otherObject.prevPageCursor)
        .isEquals();
  }
}
