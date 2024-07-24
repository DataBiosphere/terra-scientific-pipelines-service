package bio.terra.pipelines.common.utils.pagination;

/**
 * Class used to structure the response from a paginated call to the database
 *
 * @param content - records returned from database
 * @param previousPageCursor - if exists, the page token for the next page
 * @param nextPageCursor - if exists, the page token for the previous page
 * @param <T> - Db entity queried over
 */
public record PageResponse<T>(T content, String previousPageCursor, String nextPageCursor) {}
