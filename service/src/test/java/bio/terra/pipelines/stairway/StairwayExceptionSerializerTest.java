package bio.terra.pipelines.stairway;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.dependencies.common.DependencyNotAvailableException;
import bio.terra.pipelines.dependencies.stairway.StairwayExceptionSerializer;
import bio.terra.pipelines.dependencies.stairway.exception.ExceptionSerializerException;
import bio.terra.pipelines.testutils.BaseTest;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.apache.commons.lang3.StringUtils;
import org.junit.jupiter.api.Test;

class StairwayExceptionSerializerTest extends BaseTest {
  private final ObjectMapper objectMapper = new ObjectMapper();

  StairwayExceptionSerializer stairwayExceptionSerializer =
      new StairwayExceptionSerializer(objectMapper);

  @Test
  void serialize_nullException() {
    assertEquals(StringUtils.EMPTY, stairwayExceptionSerializer.serialize(null));
  }

  @Test
  void serialize_nonRuntimeException() {
    Exception exception = new Exception("test exception");
    String serializedException = stairwayExceptionSerializer.serialize(exception);
    String expected =
        "{\"className\":\"bio.terra.pipelines.dependencies.stairway.exception.StairwayJobResponseException\",\"message\":\"test exception\",\"errorDetails\":[],\"errorCode\":500,\"apiErrorReportException\":true}";
    assertEquals(expected, serializedException);
  }

  @Test
  void serialize_errorReportException() {
    RuntimeException exception =
        new DependencyNotAvailableException("test dependency name", "test context");
    String serializedException = stairwayExceptionSerializer.serialize(exception);
    String expected =
        "{\"className\":\"bio.terra.pipelines.dependencies.common.DependencyNotAvailableException\",\"message\":\"Dependency not available: test dependency name. test context\",\"errorDetails\":[],\"errorCode\":500,\"apiErrorReportException\":true}";
    assertEquals(expected, serializedException);
  }

  @Test
  void serialize_nonErrorReportException() {
    RuntimeException exception = new RuntimeException("test exception");
    String serializedException = stairwayExceptionSerializer.serialize(exception);
    String expected =
        "{\"className\":\"java.lang.RuntimeException\",\"message\":\"test exception\",\"errorDetails\":null,\"errorCode\":0,\"apiErrorReportException\":false}";
    assertEquals(expected, serializedException);
  }

  @Test
  void deserialize_nullException() {
    assertEquals(null, stairwayExceptionSerializer.deserialize(null));
  }

  @Test
  void deserialize_badInput() {
    String badInput = "bad input";

    Exception deserializedException = stairwayExceptionSerializer.deserialize(badInput);

    assertEquals(ExceptionSerializerException.class, deserializedException.getClass());
    assertTrue(
        deserializedException.getMessage().contains("Failed to deserialize exception data:"),
        "Exception message should contain 'Failed to deserialize exception data:'");
  }

  @Test
  void deserialize_classNotFound() {
    String serializedException =
        "{\"className\":\"unrecognized class name\",\"message\":\"test exception\",\"errorDetails\":[],\"errorCode\":500,\"apiErrorReportException\":true}";

    Exception deserializedException = stairwayExceptionSerializer.deserialize(serializedException);

    assertTrue(
        deserializedException
            .getMessage()
            .contains("Exception class not found: unrecognized class name"),
        "Exception message should contain 'Exception class not found: [unrecognized class name]'");
  }

  @Test
  void deserialize_failToConstruct() {
    // this is marked as an ApiErrorReportException, but no matching constructor can be found
    String serializedException =
        "{\"className\":\"bio.terra.pipelines.dependencies.common.DependencyNotAvailableException\",\"message\":\"test message\",\"errorDetails\":[],\"errorCode\":500,\"apiErrorReportException\":true}";

    Exception deserializedException = stairwayExceptionSerializer.deserialize(serializedException);

    assertEquals(ExceptionSerializerException.class, deserializedException.getClass());
    assertEquals(
        "Failed to construct exception: bio.terra.pipelines.dependencies.common.DependencyNotAvailableException; Exception message: test message",
        deserializedException.getMessage());
  }
}
