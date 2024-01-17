package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.common.exception.MissingRequiredFieldsException;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.dependencies.stairway.exception.InvalidResultStateException;
import bio.terra.pipelines.generated.model.ApiErrorReport;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import com.fasterxml.jackson.core.type.TypeReference;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.http.HttpStatus;

class FlightUtilsTest extends BaseEmbeddedDbTest {

  @Mock private FlightContext flightContext;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    var workingMap = new FlightMap();

    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void setErrorResponse_success() {
    String message = "message";
    FlightUtils.setErrorResponse(flightContext, message, HttpStatus.I_AM_A_TEAPOT);

    FlightMap workingMap = flightContext.getWorkingMap();
    ApiErrorReport response =
        workingMap.get(JobMapKeys.RESPONSE.getKeyName(), ApiErrorReport.class);

    assertNotNull(response);
    assertEquals(message, response.getMessage());
    assertEquals(
        HttpStatus.I_AM_A_TEAPOT,
        workingMap.get(JobMapKeys.STATUS_CODE.getKeyName(), HttpStatus.class));
  }

  @Test
  void getInputParameterOrWorkingValue_fromInputs() {
    // put the kvp in the input parameters
    String key = "key";
    String value = "value";
    FlightMap inputParameters = flightContext.getInputParameters();
    inputParameters.put(key, value);

    assertEquals(
        value, FlightUtils.getInputParameterOrWorkingValue(flightContext, key, key, String.class));
  }

  @Test
  void getInputParameterOrWorkingValue_fromWorkingMap() {
    // put the kvp in the working map but not in the input parameters
    String key = "key";
    String value = "value";
    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(key, value);

    assertEquals(
        value, FlightUtils.getInputParameterOrWorkingValue(flightContext, key, key, String.class));
  }

  @Test
  void getInputParametersOrWorkingValue_null() {
    // put the kvp in the working map but not in the input parameters
    String key = "key";
    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(key, null);

    assertNull(FlightUtils.getInputParameterOrWorkingValue(flightContext, key, key, String.class));
  }

  @Test
  void validateRequiredEntries_success() {
    String requiredKey1 = "requiredKey1";
    String requiredKey2 = "requiredKey2";
    FlightMap flightMap = new FlightMap();

    flightMap.put(requiredKey1, "value1");
    flightMap.put(requiredKey2, "value2");

    assertDoesNotThrow(
        () -> FlightUtils.validateRequiredEntries(flightMap, requiredKey1, requiredKey2));
  }

  @Test
  void validateRequiredEntries_missingRequiredKey() {
    String requiredKey1 = "requiredKey1";
    String requiredKey2 = "requiredKey2";
    FlightMap flightMap = new FlightMap();

    flightMap.put(requiredKey1, "value1");

    assertThrows(
        MissingRequiredFieldsException.class,
        () -> FlightUtils.validateRequiredEntries(flightMap, requiredKey1, requiredKey2));
  }

  @Test
  void getResultMapRequired_success() {
    FlightState flightState = new FlightState();
    FlightMap resultMap = new FlightMap();
    flightState.setResultMap(resultMap);

    FlightMap resultFromFlightState = FlightUtils.getResultMapRequired(flightState);
    assertEquals(resultMap, resultFromFlightState);
  }

  @Test
  void getResultMapRequired_missingResultMap() {
    FlightState flightState = new FlightState();

    assertThrows(
        MissingRequiredFieldsException.class, () -> FlightUtils.getResultMapRequired(flightState));
  }

  @Test
  void getFlightErrorMessage_withMessage() {
    FlightState flightState = new FlightState();
    String message = "message";

    // the exact exception type doesn't matter
    InvalidResultStateException exception = new InvalidResultStateException(message);
    flightState.setException(exception);

    assertEquals(message, FlightUtils.getFlightErrorMessage(flightState));
  }

  @Test
  void getFlightErrorMessage_noMessage() {
    FlightState flightState = new FlightState();

    // the exact exception type doesn't matter, but had to find one that accepts no message
    IndexOutOfBoundsException exception = new IndexOutOfBoundsException();
    flightState.setException(exception);

    assertEquals(
        "Exception: java.lang.IndexOutOfBoundsException",
        FlightUtils.getFlightErrorMessage(flightState));
  }

  @Test
  void getRequired_class_success() {
    String key = "key";
    String value = "value";
    FlightMap flightMap = new FlightMap();
    flightMap.put(key, value);

    assertEquals(value, FlightUtils.getRequired(flightMap, key, String.class));
  }

  @Test
  void getRequired_class_fail() {
    FlightMap flightMap = new FlightMap();

    assertThrows(
        MissingRequiredFieldsException.class,
        () -> FlightUtils.getRequired(flightMap, "key", String.class));
  }

  @Test
  void getRequired_typeRef_success() {
    String key = "key";
    String value = "value";
    FlightMap flightMap = new FlightMap();
    flightMap.put(key, value);

    assertEquals(value, FlightUtils.getRequired(flightMap, key, new TypeReference<>() {}));
  }

  @Test
  void getRequired_typeRef_fail() {
    FlightMap flightMap = new FlightMap();

    TypeReference<Object> typeReference = new TypeReference<>() {};

    assertThrows(
        MissingRequiredFieldsException.class,
        () -> FlightUtils.getRequired(flightMap, "key", typeReference));
  }

  @Test
  void flightComplete_success_isComplete() {
    FlightState flightState = new FlightState();
    flightState.setFlightStatus(FlightStatus.SUCCESS);

    // flightComplete returns True if flight is SUCCESS, ERROR, or FATAL
    assertTrue(FlightUtils.flightComplete(flightState));
  }

  @Test
  void flightComplete_error_isComplete() {
    FlightState flightState = new FlightState();
    flightState.setFlightStatus(FlightStatus.ERROR);

    // flightComplete returns True if flight is SUCCESS, ERROR, or FATAL
    assertTrue(FlightUtils.flightComplete(flightState));
  }

  @Test
  void flightComplete_fatal_isComplete() {
    FlightState flightState = new FlightState();
    flightState.setFlightStatus(FlightStatus.FATAL);

    // flightComplete returns True if flight is SUCCESS, ERROR, or FATAL
    assertTrue(FlightUtils.flightComplete(flightState));
  }

  @Test
  void flightComplete_running_isNotComplete() {
    FlightState flightState = new FlightState();
    flightState.setFlightStatus(FlightStatus.RUNNING);

    // flightComplete returns True if flight is SUCCESS, ERROR, or FATAL
    assertFalse(FlightUtils.flightComplete(flightState));
  }
}
