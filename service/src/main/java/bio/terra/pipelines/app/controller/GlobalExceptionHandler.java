package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.ErrorReportException;
import bio.terra.pipelines.generated.model.ApiErrorReport;
import io.sentry.Sentry;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import javax.validation.constraints.NotNull;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.http.converter.HttpMessageNotReadableException;
import org.springframework.retry.backoff.BackOffInterruptedException;
import org.springframework.web.HttpRequestMethodNotSupportedException;
import org.springframework.web.bind.MethodArgumentNotValidException;
import org.springframework.web.bind.annotation.ExceptionHandler;
import org.springframework.web.bind.annotation.RestControllerAdvice;
import org.springframework.web.method.annotation.MethodArgumentTypeMismatchException;
import org.springframework.web.servlet.NoHandlerFoundException;

// This module provides a top-level exception handler for controllers.
// All exceptions that rise through the controllers are caught in this handler.
// It converts the exceptions into standard ApiErrorReport responses.

@RestControllerAdvice
public class GlobalExceptionHandler {
  private final Logger logger = LoggerFactory.getLogger(GlobalExceptionHandler.class);

  // -- Error Report - one of our exceptions --
  @ExceptionHandler(ErrorReportException.class)
  public ResponseEntity<ApiErrorReport> errorReportHandler(ErrorReportException ex) {
    return buildApiErrorReport(ex, ex.getStatusCode(), ex.getCauses());
  }

  // -- validation exceptions - we don't control the exception raised
  @ExceptionHandler({
    MethodArgumentTypeMismatchException.class,
    HttpRequestMethodNotSupportedException.class,
    IllegalArgumentException.class,
    NoHandlerFoundException.class
  })
  public ResponseEntity<ApiErrorReport> validationExceptionHandler(Exception ex) {
    logger.error("Global exception handler: catch stack", ex);
    // For security reasons, we generally don't want to include the user's invalid (and potentially
    // malicious) input in the error response, which also means we don't include the full exception.
    // Instead, we return a generic error message about input validation.
    String validationErrorMessage =
        "Request could not be parsed or was invalid: "
            + ex.getClass().getSimpleName()
            + ". Ensure that all types are correct and that enums have valid values.";
    ApiErrorReport errorReport =
        new ApiErrorReport()
            .message(validationErrorMessage)
            .statusCode(HttpStatus.BAD_REQUEST.value());
    return new ResponseEntity<>(errorReport, HttpStatus.BAD_REQUEST);
  }

  @ExceptionHandler({MethodArgumentNotValidException.class})
  public ResponseEntity<ApiErrorReport> argumentNotValidExceptionHandler(
      MethodArgumentNotValidException ex) {
    logger.debug("MethodArgumentNotValid exception caught by global exception handler", ex);
    // For security reasons, we generally don't want to include the user's invalid (and potentially
    // malicious) input in the error response, which also means we don't include the full exception.
    // Instead, we return a generic error message about input validation.
    List<String> errors = new ArrayList<>();
    ex.getBindingResult()
        .getFieldErrors()
        .forEach(
            error -> {
              String fieldName = error.getField();
              String errorMessage = error.getDefaultMessage();
              errors.add(fieldName + " " + errorMessage);
            });
    Collections.sort(errors); // sort alphabetically to make testing easier
    String validationErrorMessage =
        "Request could not be parsed or was invalid: " + String.join("; ", errors);
    ApiErrorReport errorReport =
        new ApiErrorReport()
            .message(validationErrorMessage)
            .statusCode(HttpStatus.BAD_REQUEST.value());
    return new ResponseEntity<>(errorReport, HttpStatus.BAD_REQUEST);
  }

  @ExceptionHandler({HttpMessageNotReadableException.class})
  public ResponseEntity<ApiErrorReport> httpMessageNotReadableExceptionHandler(
      HttpMessageNotReadableException ex) {
    logger.debug(
        "HttpMessageNotReadableException exception caught by global exception handler", ex);

    // Extract the top-level error message without the nested exceptions
    String message = ex.getMessage();
    String validationErrorMessage = message;
    final int tailIndex = StringUtils.indexOf(message, "; nested exception is");
    if (tailIndex != -1) {
      validationErrorMessage = StringUtils.left(message, tailIndex);
    }

    ApiErrorReport errorReport =
        new ApiErrorReport()
            .message(validationErrorMessage)
            .statusCode(HttpStatus.BAD_REQUEST.value());
    return new ResponseEntity<>(errorReport, HttpStatus.BAD_REQUEST);
  }

  // Note we don't use @retryable yet, but we anticipate using it later so leaving it for now
  // Exception thrown by Spring Retry code when interrupted.
  @ExceptionHandler({BackOffInterruptedException.class})
  public ResponseEntity<ApiErrorReport> retryBackoffExceptionHandler(
      BackOffInterruptedException ex) {
    String errorMessage =
        "Unexpected interrupt while retrying database logic. This may succeed on a retry. "
            + ex.getMessage();
    ApiErrorReport errorReport =
        new ApiErrorReport()
            .message(errorMessage)
            .statusCode(HttpStatus.INTERNAL_SERVER_ERROR.value());
    return new ResponseEntity<>(errorReport, HttpStatus.INTERNAL_SERVER_ERROR);
  }

  // -- catchall - log so we can understand what we have missed in the handlers above
  @ExceptionHandler(Exception.class)
  public ResponseEntity<ApiErrorReport> catchallHandler(Exception ex) {
    logger.error("Exception caught by catchall handler", ex);
    return buildApiErrorReport(ex, HttpStatus.INTERNAL_SERVER_ERROR, null);
  }

  private ResponseEntity<ApiErrorReport> buildApiErrorReport(
      @NotNull Throwable ex, HttpStatus statusCode, List<String> causes) {
    // only logging 5** errors to sentry, should we log more/every error?  Other services seem to do 5** only
    if (statusCode.is5xxServerError()) {
      Sentry.captureException(ex);
    }
    StringBuilder combinedCauseString = new StringBuilder();
    for (Throwable cause = ex; cause != null; cause = cause.getCause()) {
      combinedCauseString.append("cause: ").append(cause).append(", ");
    }
    logger.error("Global exception handler: {}", combinedCauseString, ex);

    ApiErrorReport errorReport =
        new ApiErrorReport().message(ex.getMessage()).statusCode(statusCode.value()).causes(causes);
    return new ResponseEntity<>(errorReport, statusCode);
  }
}
