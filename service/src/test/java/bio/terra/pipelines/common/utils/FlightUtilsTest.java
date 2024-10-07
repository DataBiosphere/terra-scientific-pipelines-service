package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;

import bio.terra.common.exception.MissingRequiredFieldException;
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
  void setErrorResponseSuccess() {
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
  void getInputParameterOrWorkingValueFromInputs() {
    // put the kvp in the input parameters
    String key = "key";
    String value = "value";
    FlightMap inputParameters = flightContext.getInputParameters();
    inputParameters.put(key, value);

    assertEquals(
        value, FlightUtils.getInputParameterOrWorkingValue(flightContext, key, key, String.class));
  }

  @Test
  void getInputParameterOrWorkingValueFromWorkingMap() {
    // put the kvp in the working map but not in the input parameters
    String key = "key";
    String value = "value";
    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(key, value);

    assertEquals(
        value, FlightUtils.getInputParameterOrWorkingValue(flightContext, key, key, String.class));
  }

  @Test
  void getInputParametersOrWorkingValueNull() {
    // put the kvp in the working map but not in the input parameters
    String key = "key";
    FlightMap workingMap = flightContext.getWorkingMap();
    workingMap.put(key, null);

    assertNull(FlightUtils.getInputParameterOrWorkingValue(flightContext, key, key, String.class));
  }

  @Test
  void validateRequiredEntriesSuccess() {
    String requiredKey1 = "requiredKey1";
    String requiredKey2 = "requiredKey2";
    FlightMap flightMap = new FlightMap();

    flightMap.put(requiredKey1, "value1");
    flightMap.put(requiredKey2, "value2");

    assertDoesNotThrow(
        () -> FlightUtils.validateRequiredEntries(flightMap, requiredKey1, requiredKey2));
  }

  @Test
  void validateRequiredEntriesMissingRequiredKey() {
    String requiredKey1 = "requiredKey1";
    String requiredKey2 = "requiredKey2";
    FlightMap flightMap = new FlightMap();

    flightMap.put(requiredKey1, "value1");

    assertThrows(
        MissingRequiredFieldException.class,
        () -> FlightUtils.validateRequiredEntries(flightMap, requiredKey1, requiredKey2));
  }

  @Test
  void getResultMapRequiredSuccess() {
    FlightState flightState = new FlightState();
    FlightMap resultMap = new FlightMap();
    flightState.setResultMap(resultMap);

    FlightMap resultFromFlightState = FlightUtils.getResultMapRequired(flightState);
    assertEquals(resultMap, resultFromFlightState);
  }

  @Test
  void getResultMapRequiredMissingResultMap() {
    FlightState flightState = new FlightState();

    assertThrows(
        MissingRequiredFieldException.class, () -> FlightUtils.getResultMapRequired(flightState));
  }

  @Test
  void getFlightErrorMessageWithMessage() {
    FlightState flightState = new FlightState();
    String message = "message";

    // the exact exception type doesn't matter
    InvalidResultStateException exception = new InvalidResultStateException(message);
    flightState.setException(exception);

    assertEquals(message, FlightUtils.getFlightErrorMessage(flightState));
  }

  @Test
  void getFlightErrorMessageNoMessage() {
    FlightState flightState = new FlightState();

    // the exact exception type doesn't matter, but had to find one that accepts no message
    IndexOutOfBoundsException exception = new IndexOutOfBoundsException();
    flightState.setException(exception);

    assertEquals(
        "Exception: java.lang.IndexOutOfBoundsException",
        FlightUtils.getFlightErrorMessage(flightState));
  }

  @Test
  void getRequiredClassSuccess() {
    String key = "key";
    String value = "value";
    FlightMap flightMap = new FlightMap();
    flightMap.put(key, value);

    assertEquals(value, FlightUtils.getRequired(flightMap, key, String.class));
  }

  @Test
  void getRequiredClassFail() {
    FlightMap flightMap = new FlightMap();

    assertThrows(
        MissingRequiredFieldException.class,
        () -> FlightUtils.getRequired(flightMap, "key", String.class));
  }

  @Test
  void getRequiredTypeRefSuccess() {
    String key = "key";
    String value = "value";
    FlightMap flightMap = new FlightMap();
    flightMap.put(key, value);

    assertEquals(value, FlightUtils.getRequired(flightMap, key, new TypeReference<>() {}));
  }

  @Test
  void getRequiredTypeRefFail() {
    FlightMap flightMap = new FlightMap();

    TypeReference<Object> typeReference = new TypeReference<>() {};

    assertThrows(
        MissingRequiredFieldException.class,
        () -> FlightUtils.getRequired(flightMap, "key", typeReference));
  }

  @Test
  void flightCompleteSuccessIsComplete() {
    FlightState flightState = new FlightState();
    flightState.setFlightStatus(FlightStatus.SUCCESS);

    // flightComplete returns True if flight is SUCCESS, ERROR, or FATAL
    assertTrue(FlightUtils.flightComplete(flightState));
  }

  @Test
  void flightCompleteErrorIsComplete() {
    FlightState flightState = new FlightState();
    flightState.setFlightStatus(FlightStatus.ERROR);

    // flightComplete returns True if flight is SUCCESS, ERROR, or FATAL
    assertTrue(FlightUtils.flightComplete(flightState));
  }

  @Test
  void flightCompleteFatalIsComplete() {
    FlightState flightState = new FlightState();
    flightState.setFlightStatus(FlightStatus.FATAL);

    // flightComplete returns True if flight is SUCCESS, ERROR, or FATAL
    assertTrue(FlightUtils.flightComplete(flightState));
  }

  @Test
  void flightCompleteRunningIsNotComplete() {
    FlightState flightState = new FlightState();
    flightState.setFlightStatus(FlightStatus.RUNNING);

    // flightComplete returns True if flight is SUCCESS, ERROR, or FATAL
    assertFalse(FlightUtils.flightComplete(flightState));
  }

  @Test
  void inputParametersContainTrue() {
    // key is present, value is true
    FlightMap inputParameters1 = flightContext.getInputParameters();
    inputParameters1.put("key", true);
    assertTrue(FlightUtils.inputParametersContainTrue(inputParameters1, "key"));

    // key is present, value is false
    FlightMap inputParameters2 = flightContext.getInputParameters();
    inputParameters2.put("key", false);
    assertFalse(FlightUtils.inputParametersContainTrue(inputParameters2, "key"));

    // key is not present
    FlightMap inputParameters3 = flightContext.getInputParameters();
    assertFalse(FlightUtils.inputParametersContainTrue(inputParameters3, "key"));
  }
}
