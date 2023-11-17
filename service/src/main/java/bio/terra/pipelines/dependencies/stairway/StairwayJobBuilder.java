package bio.terra.pipelines.dependencies.stairway;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.MissingRequiredFieldException;
import bio.terra.common.stairway.StairwayComponent;
import bio.terra.pipelines.common.utils.MdcHook;
import bio.terra.stairway.Flight;
import bio.terra.stairway.FlightMap;
import com.fasterxml.jackson.core.type.TypeReference;
import io.opencensus.contrib.spring.aop.Traced;
import javax.annotation.Nullable;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class StairwayJobBuilder {
  private final StairwayJobService stairwayJobService;
  private final StairwayComponent stairwayComponent;
  private final MdcHook mdcHook;
  private final FlightMap jobParameterMap;
  private Class<? extends Flight> flightClass;
  private final Logger logger = LoggerFactory.getLogger(StairwayJobBuilder.class);
  @Nullable private String jobId;
  @Nullable private String description;

  public StairwayJobBuilder(
      StairwayJobService stairwayJobService, StairwayComponent stairwayComponent, MdcHook mdcHook) {
    this.stairwayJobService = stairwayJobService;
    this.stairwayComponent = stairwayComponent;
    this.mdcHook = mdcHook;
    this.jobParameterMap = new FlightMap();
  }

  public enum JobMapKeys {
    // parameters for all flight types
    DESCRIPTION("description"),
    RESPONSE("response"),
    STATUS_CODE("status_code"),
    RESULT_PATH("resultPath"),

    // parameter for the job
    FLIGHT_CLASS("flight_class");

    private final String keyName;

    JobMapKeys(String keyName) {
      this.keyName = keyName;
    }

    public String getKeyName() {
      return keyName;
    }

    public static boolean isRequiredKey(String keyName) {
      return keyName.equals(JobMapKeys.DESCRIPTION.getKeyName());
    }
  }

  public StairwayJobBuilder flightClass(Class<? extends Flight> flightClass) {
    this.flightClass = flightClass;
    return this;
  }

  public StairwayJobBuilder stairwayJobId(@Nullable String stairwayJobId) {
    // If clients provide a non-null job ID, it cannot be whitespace-only
    if (StringUtils.isWhitespace(stairwayJobId)) {
      throw new BadRequestException("stairwayJobId cannot be whitespace-only.");
    }
    this.jobId = stairwayJobId;
    return this;
  }

  public StairwayJobBuilder description(@Nullable String description) {
    this.description = description;
    return this;
  }

  public StairwayJobBuilder addParameter(String keyName, @Nullable Object val) {
    if (StringUtils.isBlank(keyName)) {
      throw new BadRequestException("Parameter name cannot be null or blanks.");
    }
    // note that this call overwrites a parameter if it already exists
    jobParameterMap.put(keyName, val);
    return this;
  }

  /**
   * Submit a job to stairway and return the jobID immediately.
   *
   * @return jobID of submitted flight
   */
  public String submit() {
    populateInputParams();
    return stairwayJobService.submit(flightClass, jobParameterMap, jobId);
  }

  /**
   * Submit a job to stairway, wait until it's complete, and return the job result.
   *
   * @param resultClass Class of the job's result
   * @return Result of the finished job.
   */
  @Traced
  public <T> T submitAndWait(Class<T> resultClass) {
    populateInputParams();
    return stairwayJobService.submitAndWait(
        flightClass, jobParameterMap, resultClass, /*typeReference=*/ null, jobId);
  }

  /**
   * Submit a job to stairway, wait until it's complete, and return the job result.
   *
   * @param typeReference Class of the job's result
   * @return Result of the finished job.
   */
  @Traced
  public <T> T submitAndWait(TypeReference<T> typeReference) {
    populateInputParams();
    return stairwayJobService.submitAndWait(
        flightClass, jobParameterMap, /*resultClass=*/ null, typeReference, jobId);
  }

  /**
   * Submit a job to stairway, wait until it's complete, and return the job result.
   *
   * @return Result of the finished job.
   */
  @Traced
  public <T> T submitAndWait() {
    populateInputParams();
    return stairwayJobService.submitAndWait(
        flightClass, jobParameterMap, /*resultClass=*/ null, /*typeReference=*/ null, jobId);
  }

  // Check the inputs, supply defaults and finalize the input parameter map
  private void populateInputParams() {
    if (flightClass == null) {
      throw new MissingRequiredFieldException("Missing flight class: flightClass");
    }

    // Default to a generated job id
    if (jobId == null) {
      jobId = stairwayComponent.get().createFlightId();
    }

    // Always add the MDC logging and tracing span parameters for the mdc hook
    addParameter(MdcHook.MDC_FLIGHT_MAP_KEY, mdcHook.getSerializedCurrentContext());
    //        addParameter(
    //                MonitoringHook.SUBMISSION_SPAN_CONTEXT_MAP_KEY,
    //                MonitoringHook.serializeCurrentTracingContext());

    // Convert any other members that were set into parameters. However, if they were
    // explicitly added with addParameter during construction, we do not overwrite them.
    if (shouldInsert(JobMapKeys.DESCRIPTION, description)) {
      addParameter(JobMapKeys.DESCRIPTION.getKeyName(), description);
    }
  }

  private boolean shouldInsert(String mapKey, @Nullable Object value) {
    return (value != null && !jobParameterMap.containsKey(mapKey));
  }

  private boolean shouldInsert(JobMapKeys mapKey, @Nullable Object value) {
    return (value != null && !jobParameterMap.containsKey(mapKey.getKeyName()));
  }
}
