package bio.terra.pipelines.dependencies.stairway;

import bio.terra.common.db.DataSourceInitializer;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.common.logging.LoggingUtils;
import bio.terra.common.stairway.MonitoringHook;
import bio.terra.common.stairway.StairwayComponent;
import bio.terra.pipelines.app.configuration.internal.StairwayDatabaseConfiguration;
import bio.terra.pipelines.app.configuration.internal.StairwayJobConfiguration;
import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.MdcHook;
import bio.terra.pipelines.dependencies.stairway.exception.*;
import bio.terra.stairway.*;
import bio.terra.stairway.exception.DuplicateFlightIdException;
import bio.terra.stairway.exception.FlightNotFoundException;
import bio.terra.stairway.exception.StairwayException;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.introspect.AnnotatedMember;
import com.fasterxml.jackson.databind.introspect.JacksonAnnotationIntrospector;
import com.google.common.annotations.VisibleForTesting;
import io.opencensus.contrib.spring.aop.Traced;
import java.time.Duration;
import java.util.Optional;
import javax.annotation.Nullable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class StairwayJobService {
  private static final int METADATA_ROW_WAIT_SECONDS = 1;
  private static final Duration METADATA_ROW_MAX_WAIT_TIME = Duration.ofSeconds(28);

  private final StairwayJobConfiguration stairwayJobConfig;
  private final StairwayDatabaseConfiguration stairwayDatabaseConfig;
  private final MdcHook mdcHook;
  private final StairwayComponent stairwayComponent;
  private final FlightBeanBag flightBeanBag;
  private final Logger logger = LoggerFactory.getLogger(StairwayJobService.class);
  private final ObjectMapper objectMapper;
  private FlightDebugInfo flightDebugInfo;

  @Autowired
  public StairwayJobService(
      StairwayJobConfiguration stairwayJobConfig,
      StairwayDatabaseConfiguration stairwayDatabaseConfig,
      MdcHook mdcHook,
      StairwayComponent stairwayComponent,
      FlightBeanBag flightBeanBag,
      ObjectMapper objectMapper) {
    this.stairwayJobConfig = stairwayJobConfig;
    this.stairwayDatabaseConfig = stairwayDatabaseConfig;
    this.mdcHook = mdcHook;
    this.stairwayComponent = stairwayComponent;
    this.flightBeanBag = flightBeanBag;
    this.objectMapper = objectMapper;
  }

  // Fully fluent style of JobBuilder
  public StairwayJobBuilder newJob() {
    return new StairwayJobBuilder(this, stairwayComponent, mdcHook);
  }

  // submit a new job to stairway
  // protected method intended to be called only from JobBuilder
  protected String submit(
      Class<? extends Flight> flightClass, FlightMap parameterMap, String jobId) {
    try {
      stairwayComponent
          .get()
          .submitWithDebugInfo(
              jobId, flightClass, parameterMap, /* shouldQueue= */ false, flightDebugInfo);
    } catch (DuplicateFlightIdException ex) {
      // DuplicateFlightIdException is a more specific StairwayException, and so needs to
      // be checked separately. Allowing duplicate FlightIds is useful for ensuring idempotent
      // behavior of flights.
      logger.warn("Received duplicate job ID: {}", jobId);
      throw new DuplicateStairwayJobIdException(
          String.format("Received duplicate jobId %s", jobId), ex);
    } catch (StairwayException | InterruptedException stairwayEx) {
      throw new InternalStairwayException(stairwayEx);
    }
    return jobId;
  }

  // Submit a new job to stairway, wait for it to finish, then return the result.
  // This will throw any exception raised by the flight.
  // protected method intended to be called only from JobBuilder
  protected <T> T submitAndWait(
      Class<? extends Flight> flightClass,
      FlightMap parameterMap,
      Class<T> resultClass,
      TypeReference<T> typeReference,
      String jobId) {
    submit(flightClass, parameterMap, jobId);
    waitForJob(jobId);

    JobResultOrException<T> resultOrException =
        retrieveJobResult(jobId, resultClass, typeReference);
    if (resultOrException.getException() != null) {
      throw resultOrException.getException();
    }
    return resultOrException.getResult();
  }

  /**
   * This method is called from StartupInitializer as part of the sequence of migrating databases
   * and recovering any jobs; i.e., Stairway flights. It is moved here so that JobService
   * encapsulates all of the Stairway interaction.
   */
  public void initialize() {
    configureMapper();
    stairwayComponent.initialize(
        stairwayComponent
            .newStairwayOptionsBuilder()
            .dataSource(DataSourceInitializer.initializeDataSource(stairwayDatabaseConfig))
            .context(flightBeanBag)
            .addHook(mdcHook)
            .addHook(new MonitoringHook())
            .exceptionSerializer(new StairwayExceptionSerializer(objectMapper)));
  }

  /**
   * This is currently a hack, because Stairway does not provide a way to pass in the mapper. It
   * does expose its own mapper for testing, so we use that public API to add the introspector that
   * we need.
   *
   * <p>TODO: PF-2505 When that Stairway feature is done we should create and set our own object
   * mapper in Stairway.
   */
  @VisibleForTesting
  public static void configureMapper() {
    StairwayMapper.getObjectMapper().setAnnotationIntrospector(new IgnoreInheritedIntrospector());
  }

  /**
   * Jackson does not see @JsonIgnore annotations from super classes. That means any getter in a
   * super class gets serialized by default. We do not want that behavior, so we add this
   * introspector to the ObjectMapper to force ignore everything from the resource super classes.
   *
   * <p>We do not need to ignore ReferencedResource class because it does not add any fields to
   * those coming from WsmResource class.
   */
  private static class IgnoreInheritedIntrospector extends JacksonAnnotationIntrospector {
    @Override
    public boolean hasIgnoreMarker(AnnotatedMember m) {
      // TODO see if we need to adapt this for TSPS
      //      boolean ignore =
      //              (m.getDeclaringClass() == WsmResource.class)
      //                      || (m.getDeclaringClass() == ControlledResource.class);
      //      return ignore || super.hasIgnoreMarker(m);
      return super.hasIgnoreMarker(m);
    }
  }

  public void waitForJob(String jobId) {
    try {
      FlightUtils.waitForJobFlightCompletion(stairwayComponent.get(), jobId);
    } catch (Exception ex) {
      throw new InternalServerErrorException(ex);
    }
  }

  /** Retrieves Job Result specifying the result class type. */
  @Traced
  public <T> JobResultOrException<T> retrieveJobResult(String jobId, Class<T> resultClass) {
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
   *       we store StairwayJobResponseException, passing in the Throwable to the exception. In the
   *       global exception handler, we retrieve the Throwable and use the error text from that in
   *       the error model.
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
      String jobId, @Nullable Class<T> resultClass, @Nullable TypeReference<T> typeReference) {
    try {
      FlightState flightState = stairwayComponent.get().getFlightState(jobId);
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
                .result(resultMap.get(StairwayJobMapKeys.RESPONSE.getKeyName(), resultClass));
          }
          if (typeReference != null) {
            return new JobResultOrException<T>()
                .result(resultMap.get(StairwayJobMapKeys.RESPONSE.getKeyName(), typeReference));
          }
          return new JobResultOrException<T>()
              .result(resultMap.get(StairwayJobMapKeys.RESPONSE.getKeyName(), (Class<T>) null));
        case RUNNING:
          throw new StairwayJobNotCompleteException(
              "Attempt to retrieve job result before job is complete; job id: "
                  + flightState.getFlightId());
        default:
          throw new InvalidResultStateException("Impossible case reached");
      }
    } catch (FlightNotFoundException flightNotFoundException) {
      throw new StairwayJobNotFoundException(
          "The flight " + jobId + " was not found", flightNotFoundException);
    } catch (StairwayException | InterruptedException stairwayEx) {
      throw new InternalStairwayException(stairwayEx);
    }
  }

  // TODO do we want this annotation / to use FindBugs?
  // @SuppressFBWarnings(value = "NM_CLASS_NOT_EXCEPTION", justification = "Non-exception by
  // design.")
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
      if (exception instanceof RuntimeException) {
        return new JobResultOrException<T>().exception((RuntimeException) exception);
      } else {
        return new JobResultOrException<T>()
            .exception(new InternalServerErrorException("wrap non-runtime exception", exception));
      }
    }
    logAlert("TSPS Stairway flight {} failed with no exception given", flightState.getFlightId());
    throw new InvalidResultStateException("Failed operation with no exception reported.");
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
