package bio.terra.pipelines.common.utils;

import bio.terra.common.exception.MissingRequiredFieldException;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.generated.model.ApiErrorReport;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.FlightState;
import bio.terra.stairway.FlightStatus;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.Optional;
import javax.annotation.Nullable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.http.HttpStatus;

/** Common methods for building flights */
public final class FlightUtils {
  private static final Logger logger = LoggerFactory.getLogger(FlightUtils.class);

  private FlightUtils() {}

  /**
   * Build an error model and set it as the response. Useful if we want to return an error
   * deliberately in a step.
   *
   * @param context flight context
   * @param message error message
   * @param responseStatus status
   */
  public static void setErrorResponse(
      FlightContext context, String message, HttpStatus responseStatus) {
    ApiErrorReport errorModel = new ApiErrorReport().message(message);
    setResponse(context, errorModel, responseStatus);
  }

  /**
   * Set the response and status code in the result map.
   *
   * @param context flight context
   * @param responseObject response object to set
   * @param responseStatus status code to set
   */
  public static void setResponse(
      FlightContext context, Object responseObject, HttpStatus responseStatus) {
    FlightMap workingMap = context.getWorkingMap();
    workingMap.put(JobMapKeys.RESPONSE.getKeyName(), responseObject);
    workingMap.put(JobMapKeys.STATUS_CODE.getKeyName(), responseStatus);
  }

  /**
   * Get a supplied input value from input parameters, or, if that's missing, a default (previous)
   * value from the working map, or null.
   *
   * @param flightContext - context object for the flight, used to get the input & working maps
   * @param inputKey - key in input parameters for the supplied (override) value
   * @param workingKey - key in the working map for the previous value
   * @param klass - class of the value, e.g. String.class
   * @param <T> - type parameter corresponding to the klass
   * @return - a value from one of the two sources, or null
   */
  public static <T> T getInputParameterOrWorkingValue(
      FlightContext flightContext, String inputKey, String workingKey, Class<T> klass) {
    return Optional.ofNullable(flightContext.getInputParameters().get(inputKey, klass))
        .orElse(flightContext.getWorkingMap().get(workingKey, klass));
  }

  /**
   * Validation function, intended to be called from the top and bottom of a doStep() method in a
   * Step. Checks a list of input (or output) string keys to ensure they have a non-null value in
   * the map. For checking the input parameters map, the call can be made from a flight constructor.
   *
   * @param flightMap - either an input parameters or working map
   * @param keys - vararg of string keys to be checked
   */
  public static void validateRequiredEntries(FlightMap flightMap, String... keys) {
    for (String key : keys) {
      if (null == flightMap.getRaw(key)) {
        throw new MissingRequiredFieldException(
            String.format("Required entry with key %s missing from flight map.", key));
      }
    }
  }

  public static FlightMap getResultMapRequired(FlightState flightState) {
    return flightState
        .getResultMap()
        .orElseThrow(
            () ->
                new MissingRequiredFieldException(
                    String.format(
                        "ResultMap is missing for flight %s", flightState.getFlightId())));
  }

  /**
   * Get the error message from a FlightState's exception object if it exists. If there is no
   * message, return a message based off the exception's class name.
   *
   * @param flightState - state of subflight
   * @return - error message or null for none
   */
  @Nullable
  public static String getFlightErrorMessage(FlightState flightState) {
    String errorMessage = flightState.getException().map(Throwable::getMessage).orElse(null);
    if (null == errorMessage && flightState.getException().isPresent()) {
      // If the exception doesn't provide a message, we can scrape the class name at least.
      errorMessage =
          flightState
              .getException()
              .map(Exception::getClass)
              .map(Class::getName)
              .map(s -> "Exception: " + s)
              .orElse(null);
    }
    return errorMessage;
  }

  /**
   * Get a value from one of the flight maps and check that it is not null. If it is null, throw.
   *
   * @param flightMap input or working map
   * @param key string key to lookup in the map
   * @param tClass class to return
   * @param <T> generic
   * @return T
   */
  public static <T> T getRequired(FlightMap flightMap, String key, Class<T> tClass) {
    var value = flightMap.get(key, tClass);
    if (value == null) {
      throw new MissingRequiredFieldException("Missing required flight map key: " + key);
    }
    return value;
  }

  /**
   * Get a value from one of the flight maps and check that it is not null. If it is null, throw.
   *
   * @param flightMap input or working map
   * @param key string key to lookup in the map
   * @param typeReference Jackson type reference
   * @param <T> generic
   * @return T
   */
  public static <T> T getRequired(FlightMap flightMap, String key, TypeReference<T> typeReference) {
    var value = flightMap.get(key, typeReference);
    if (value == null) {
      throw new MissingRequiredFieldException("Missing required flight map key: " + key);
    }
    return value;
  }

  public static boolean flightComplete(FlightState flightState) {
    logger.info(
        "Testing flight {} completion; state is {}",
        flightState.getFlightId(),
        flightState.getFlightStatus());
    return (flightState.getFlightStatus() == FlightStatus.ERROR
        || flightState.getFlightStatus() == FlightStatus.FATAL
        || flightState.getFlightStatus() == FlightStatus.SUCCESS);
  }
}
