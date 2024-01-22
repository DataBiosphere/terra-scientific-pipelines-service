package bio.terra.pipelines.dependencies.stairway;

import bio.terra.common.db.DataSourceInitializer;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.common.logging.LoggingUtils;
import bio.terra.common.stairway.MonitoringHook;
import bio.terra.common.stairway.StairwayComponent;
import bio.terra.pipelines.app.configuration.internal.StairwayDatabaseConfiguration;
import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.MdcHook;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.stairway.exception.*;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJob;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.DuplicateFlightIdException;
import bio.terra.stairway.exception.FlightNotFoundException;
import bio.terra.stairway.exception.StairwayException;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.google.common.annotations.VisibleForTesting;
import io.opencensus.contrib.spring.aop.Traced;
import java.util.*;
import javax.annotation.Nullable;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class JobService {
  private final StairwayDatabaseConfiguration stairwayDatabaseConfiguration;
  private final MdcHook mdcHook;
  private final StairwayComponent stairwayComponent;
  private final FlightBeanBag flightBeanBag;
  private final Logger logger = LoggerFactory.getLogger(JobService.class);
  private final ObjectMapper objectMapper;
  private FlightDebugInfo flightDebugInfo;

  private static final String JOB_NOT_FOUND_MSG = "The flight %s was not found";
  private static final String INTERRUPTED_MSG = "Interrupted while submitting job {}";

  @Autowired
  public JobService(
      StairwayDatabaseConfiguration stairwayDatabaseConfiguration,
      MdcHook mdcHook,
      StairwayComponent stairwayComponent,
      FlightBeanBag flightBeanBag,
      ObjectMapper objectMapper) {
    this.stairwayDatabaseConfiguration = stairwayDatabaseConfiguration;
    this.mdcHook = mdcHook;
    this.stairwayComponent = stairwayComponent;
    this.flightBeanBag = flightBeanBag;
    this.objectMapper = objectMapper;
  }

  // Fully fluent style of JobBuilder
  public JobBuilder newJob() {
    return new JobBuilder(this, mdcHook);
  }

  // submit a new job to stairway
  // protected method intended to be called only from JobBuilder
  protected UUID submit(Class<? extends Flight> flightClass, FlightMap parameterMap, UUID jobId) throws DuplicateJobIdException, InternalStairwayException {
    String jobIdString = jobId.toString();
    try {
      stairwayComponent
          .get()
          .submitWithDebugInfo(
              jobIdString, flightClass, parameterMap, /* shouldQueue= */ false, flightDebugInfo);
    } catch (DuplicateFlightIdException ex) {
      // DuplicateFlightIdException is a more specific StairwayException, and so needs to
      // be checked separately. Allowing duplicate FlightIds is useful for ensuring idempotent
      // behavior of flights.
      logger.warn("Received duplicate job ID: {}", jobIdString);
      throw new DuplicateJobIdException(
          String.format("Received duplicate jobId %s", jobIdString), ex);
    } catch (StairwayException stairwayEx) {
      throw new InternalStairwayException(stairwayEx);
    } catch (InterruptedException e) {
      logger.warn(INTERRUPTED_MSG, jobIdString);
      Thread.currentThread().interrupt();
    }
    return jobId;
  }

  /**
   * This method is called from StartupInitializer as part of the sequence of migrating databases
   * and recovering any jobs; i.e., Stairway flights. It is moved here so that JobService
   * encapsulates all of the Stairway interaction.
   */
  public void initialize() {
    stairwayComponent.initialize(
        stairwayComponent
            .newStairwayOptionsBuilder()
            .dataSource(DataSourceInitializer.initializeDataSource(stairwayDatabaseConfiguration))
            .context(flightBeanBag)
            .addHook(mdcHook)
            .addHook(new MonitoringHook())
            .exceptionSerializer(new StairwayExceptionSerializer(objectMapper)));
  }

  /** Retrieves Job Result specifying the result class type. */
  @Traced
  public <T> JobResultOrException<T> retrieveJobResult(UUID jobId, Class<T> resultClass) {
    return retrieveJobResult(jobId, resultClass, /*typeReference=*/ null);
  }

  /**
   * There are four cases to handle here:
   *
   * <ol>
   *   <li>Flight is still running. Throw an JobNotComplete exception
   *   <li>Successful flight: extract the resultMap RESPONSE as the target class.
   *   <li>Failed flight: if there is an exception, store it in the returned JobResultOrException.
   *       Note that we only store RuntimeExceptions to allow higher-level methods to throw these
   *       exceptions if they choose. Non-runtime exceptions require throw clauses on the controller
   *       methods; those are not present in the swagger-generated code, so it introduces a
   *       mismatch. Instead, in this code if the caught exception is not a runtime exception, then
   *       we store JobResponseException, passing in the Throwable to the exception. In the global
   *       exception handler, we retrieve the Throwable and use the error text from that in the
   *       error model.
   *   <li>Failed flight: no exception present. Throw an InvalidResultState exception
   * </ol>
   *
   * @param jobId to process
   * @param resultClass nullable resultClass. When not null, cast the JobResult to the given class.
   * @param typeReference nullable typeReference. When not null, cast the JobResult to generic type.
   *     When the Job does not have a result (a.k.a. null), both resultClass and typeReference are
   *     set to null.
   * @return object of the result class pulled from the result map
   */
  @Traced
  public <T> JobResultOrException<T> retrieveJobResult(
      UUID jobId, @Nullable Class<T> resultClass, @Nullable TypeReference<T> typeReference) {
    try {
      FlightState flightState = stairwayComponent.get().getFlightState(jobId.toString());
      FlightMap resultMap =
          flightState.getResultMap().orElseThrow(InvalidResultStateException::noResultMap);

      switch (flightState.getFlightStatus()) {
        case FATAL:
          logAlert("TSPS Stairway flight {} encountered dismal failure", flightState.getFlightId());
          return handleFailedFlight(flightState);
        case ERROR:
          return handleFailedFlight(flightState);
        case SUCCESS:
          if (resultClass != null) {
            return new JobResultOrException<T>()
                .result(resultMap.get(JobMapKeys.RESPONSE.getKeyName(), resultClass));
          }
          if (typeReference != null) {
            return new JobResultOrException<T>()
                .result(resultMap.get(JobMapKeys.RESPONSE.getKeyName(), typeReference));
          }
          return new JobResultOrException<T>()
              .result(resultMap.get(JobMapKeys.RESPONSE.getKeyName(), (Class<T>) null));
        case RUNNING:
          throw new JobNotCompleteException(
              "Attempt to retrieve job result before job is complete; job id: "
                  + flightState.getFlightId());
        default:
          throw new InvalidResultStateException("Impossible case reached");
      }
    } catch (FlightNotFoundException flightNotFoundException) {
      throw new JobNotFoundException(
          String.format(JOB_NOT_FOUND_MSG, jobId), flightNotFoundException);
    } catch (StairwayException stairwayEx) {
      throw new InternalStairwayException(stairwayEx);
    } catch (InterruptedException e) {
      logger.warn(INTERRUPTED_MSG, jobId);
      Thread.currentThread().interrupt();
      throw new InternalStairwayException(e);
    }
  }

  @SuppressWarnings("java:S2166") // NonExceptionNameEndsWithException by design
  public static class JobResultOrException<T> {
    private T result;
    private RuntimeException exception;

    public T getResult() {
      return result;
    }

    public JobResultOrException<T> result(T result) {
      this.result = result;
      return this;
    }

    public RuntimeException getException() {
      return exception;
    }

    public JobResultOrException<T> exception(RuntimeException exception) {
      this.exception = exception;
      return this;
    }
  }

  private <T> JobResultOrException<T> handleFailedFlight(FlightState flightState) {
    Optional<Exception> flightException = flightState.getException();
    if (flightException.isPresent()) {
      Exception exception = flightException.get();
      if (exception instanceof RuntimeException runtimeException) {
        return new JobResultOrException<T>().exception(runtimeException);
      } else {
        return new JobResultOrException<T>()
            .exception(new InternalServerErrorException("wrap non-runtime exception", exception));
      }
    }
    logAlert("TSPS Stairway flight {} failed with no exception given", flightState.getFlightId());
    throw new InvalidResultStateException("Failed operation with no exception reported.");
  }

  @VisibleForTesting
  public Stairway getStairway() {
    return stairwayComponent.get();
  }

  @Traced
  public FlightState retrieveJob(UUID jobId, String userId) {
    try {
      FlightState result = stairwayComponent.get().getFlightState(jobId.toString());
      // Note: after implementing TSPS-134, we can filter by flightId in enumerateJobs and remove
      // the following check
      if (!userId.equals(
          result.getInputParameters().get(JobMapKeys.USER_ID.getKeyName(), String.class))) {
        logger.info(
            "User {} attempted to retrieve job {} but is not the original submitter",
            userId,
            jobId);
        throw new JobUnauthorizedException(
            String.format("Caller unauthorized to access job %s", jobId));
      }
      return result;
    } catch (FlightNotFoundException flightNotFoundException) {
      throw new JobNotFoundException(
          String.format(JOB_NOT_FOUND_MSG, jobId), flightNotFoundException);
    } catch (StairwayException stairwayEx) {
      throw new InternalStairwayException(stairwayEx);
    } catch (InterruptedException e) {
      logger.warn(INTERRUPTED_MSG, jobId);
      Thread.currentThread().interrupt();
      throw new InternalStairwayException(e);
    }
  }

  /**
   * List Stairway flights submitted by a user. These inputs are translated into inputs to
   * Stairway's getFlights calls. The resulting flights are translated into enumerated jobs. The
   * jobs are ordered by submit time.
   *
   * @param userId Terra userId of the caller, to filter by
   * @param limit max number of jobs to return
   * @param pageToken optional starting place in the result set; start at beginning if missing
   * @param pipelineId optional filter by pipeline type
   * @return POJO containing the results
   */
  @Traced
  public EnumeratedJobs enumerateJobs(
      String userId, int limit, @Nullable String pageToken, @Nullable PipelinesEnum pipelineId)
      throws InternalStairwayException {
    FlightEnumeration flightEnumeration;
    try {
      FlightFilter filter = buildFlightFilter(userId, pipelineId);
      flightEnumeration = stairwayComponent.get().getFlights(pageToken, limit, filter);
    } catch (StairwayException stairwayEx) {
      throw new InternalStairwayException(stairwayEx);
    } catch (InterruptedException e) {
      logger.warn("Interrupted while enumerating jobs for user {}", userId);
      Thread.currentThread().interrupt();
      throw new InternalStairwayException(e);
    }

    List<EnumeratedJob> jobList = new ArrayList<>();
    for (FlightState state : flightEnumeration.getFlightStateList()) {
      FlightMap inputParameters = state.getInputParameters();

      String jobDescription =
          (inputParameters.containsKey(JobMapKeys.DESCRIPTION.getKeyName()))
              ? inputParameters.get(JobMapKeys.DESCRIPTION.getKeyName(), String.class)
              : StringUtils.EMPTY;

      EnumeratedJob enumeratedJob =
          new EnumeratedJob().flightState(state).jobDescription(jobDescription);
      jobList.add(enumeratedJob);
    }

    return new EnumeratedJobs()
        .pageToken(flightEnumeration.getNextPageToken())
        .totalResults(flightEnumeration.getTotalFlights())
        .results(jobList);
  }

  private FlightFilter buildFlightFilter(String userId, @Nullable PipelinesEnum pipelineId) {

    FlightFilter filter = new FlightFilter();
    // Always filter by user
    filter.addFilterInputParameter(JobMapKeys.USER_ID.getKeyName(), FlightFilterOp.EQUAL, userId);
    // Add optional filters
    Optional.ofNullable(pipelineId)
        .ifPresent(
            t ->
                filter.addFilterInputParameter(
                    JobMapKeys.PIPELINE_ID.getKeyName(), FlightFilterOp.EQUAL, t));

    return filter;
  }

  /**
   * Sets the {@link FlightDebugInfo} to manipulate future Stairway Flight submissions for testing.
   *
   * <p>This is useful for causing failures on submitted jobs. This should only be used for testing.
   */
  @VisibleForTesting
  public void setFlightDebugInfoForTest(FlightDebugInfo flightDebugInfo) {
    this.flightDebugInfo = flightDebugInfo;
  }

  private void logAlert(String msg, String flightId) {
    // Dismal and unexpected flight failures always require manual intervention,
    // so developers should be notified if they happen. The alert object is deliberately
    // not included in the error message.
    // <p>With the custom json logging configuration we have from TCL
    // (https://github.com/DataBiosphere/terra-common-lib/blob/develop/src/main/java/bio/terra/common/logging/GoogleJsonLayout.java),
    // json-like objects like alertObject in the parameter list will be included in the json
    // blob that we send to stackdriver, even if they aren't actually included in the message.
    // In this case, that means the stackdriver log will have a field terraLogBasedAlert set
    // to true which triggers alerts in the Verily deployment.
    // <p> Json objects will probably still show up in a structured way even if they're also
    // used in the log message so we could have a placeholder for it if we wanted to, but the
    // alertObject is just a map with a single value that we use to flag a message for log-based
    // alerting, so it's not particularly interesting to include in the message.
    logger.error(msg, flightId, LoggingUtils.alertObject());
  }
}
