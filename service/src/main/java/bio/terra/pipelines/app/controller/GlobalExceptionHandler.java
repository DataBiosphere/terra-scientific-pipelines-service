package bio.terra.pipelines.app.controller;

import bio.terra.common.exception.ErrorReportException;
import bio.terra.pipelines.generated.model.ApiErrorReport;
import io.sentry.Sentry;
import jakarta.annotation.Nullable;
import jakarta.validation.constraints.NotNull;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
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
import org.springframework.web.servlet.resource.NoResourceFoundException;

// This module provides a top-level exception handler for controllers.
// All exceptions that rise through the controllers are caught in this handler.
// It converts the exceptions into standard ApiErrorReport responses.

@RestControllerAdvice
public class GlobalExceptionHandler {
  private final Logger logger = LoggerFactory.getLogger(GlobalExceptionHandler.class);

  // -- Error Report - one of our exceptions --
  @ExceptionHandler({ErrorReportException.class})
  public ResponseEntity<ApiErrorReport> errorReportHandler(ErrorReportException ex) {
    return buildApiErrorReport(ex, ex.getStatusCode(), ex.getCauses(), null);
  }

  // -- validation exceptions - we don't control the exception raised
  @ExceptionHandler({
    MethodArgumentTypeMismatchException.class,
    IllegalArgumentException.class,
    NoHandlerFoundException.class
  })
  public ResponseEntity<ApiErrorReport> validationExceptionHandler(Exception ex) {
    // For security reasons, we generally don't want to include the user's invalid (and potentially
    // malicious) input in the error response, which also means we don't include the full exception.
    // Instead, we return a generic error message about input validation.
    String validationErrorMessage =
        "Request could not be parsed or was invalid: "
            + ex.getClass().getSimpleName()
            + ". Ensure that all types are correct and that enums have valid values.";
    return buildApiErrorReport(ex, HttpStatus.BAD_REQUEST, null, validationErrorMessage);
  }

  @ExceptionHandler({MethodArgumentNotValidException.class})
  public ResponseEntity<ApiErrorReport> argumentNotValidExceptionHandler(
      MethodArgumentNotValidException ex) {
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
    return buildApiErrorReport(ex, HttpStatus.BAD_REQUEST, null, validationErrorMessage);
  }

  @ExceptionHandler({HttpMessageNotReadableException.class})
  public ResponseEntity<ApiErrorReport> httpMessageNotReadableExceptionHandler(
      HttpMessageNotReadableException ex) {
    // Extract the top-level error message without the nested exceptions
    String message = ex.getMessage();
    String validationErrorMessage = message;
    final int tailIndex = StringUtils.indexOf(message, "; nested exception is");
    if (tailIndex != -1) {
      validationErrorMessage = StringUtils.left(message, tailIndex);
    }
    return buildApiErrorReport(ex, HttpStatus.BAD_REQUEST, null, validationErrorMessage);
  }

  // Note we don't use @retryable yet, but we anticipate using it later so leaving it for now
  // Exception thrown by Spring Retry code when interrupted.
  @ExceptionHandler({BackOffInterruptedException.class})
  public ResponseEntity<ApiErrorReport> retryBackoffExceptionHandler(
      BackOffInterruptedException ex) {
    String errorMessage =
        "Unexpected interrupt while retrying database logic. This may succeed on a retry. "
            + ex.getMessage();
    return buildApiErrorReport(ex, HttpStatus.INTERNAL_SERVER_ERROR, null, errorMessage);
  }

  // -- catchall - log so we can understand what we have missed in the handlers above
  @ExceptionHandler(Exception.class)
  public ResponseEntity<ApiErrorReport> catchallHandler(Exception ex) {
    logger.warn("Exception caught by catchall handler: {}", ex.getMessage());
    return buildApiErrorReport(ex, HttpStatus.INTERNAL_SERVER_ERROR, null, null);
  }

  private ResponseEntity<ApiErrorReport> buildApiErrorReport(
      @NotNull Throwable ex,
      HttpStatus statusCode,
      List<String> causes,
      @Nullable String messageForApiErrorReport) {
    // logging 4** & 5** errors to sentry
    if (statusCode.is5xxServerError() || statusCode.is4xxClientError()) {
      if (ex instanceof NoResourceFoundException
          || ex instanceof HttpRequestMethodNotSupportedException) {
        // NoResourceFoundExceptions arise from calls to nonexistent API paths and are generally
        // spam; HttpRequestMethodNotSupportedExceptions arise from calls to existing API paths with
        // unsupported methods (e.g. POST, PROPFIND) and are generally spam
        logger.warn("Not sending exception of type {} to Sentry", ex.getClass());
      } else {
        Sentry.captureException(ex);
      }
    }
    StringBuilder combinedCauseString = new StringBuilder();
    for (Throwable cause = ex; cause != null; cause = cause.getCause()) {
      combinedCauseString.append("cause: ").append(cause).append(", ");
    }
    logger.error("Global exception handler: {}", combinedCauseString, ex);

    // sanitize user facing message for 500 errors coming through controller
    if (statusCode.is5xxServerError()) {
      messageForApiErrorReport =
          "Internal server error. Please contact support if this problem persists.";
    }

    ApiErrorReport errorReport =
        new ApiErrorReport()
            .message(
                Optional.ofNullable(messageForApiErrorReport)
                    .orElse(
                        Optional.ofNullable(ex)
                            .map(Throwable::getMessage)
                            .orElse("No error message found")))
            .statusCode(statusCode.value())
            .causes(causes);
    return new ResponseEntity<>(errorReport, statusCode);
  }
}
