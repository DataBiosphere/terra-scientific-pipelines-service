package bio.terra.pipelines.dependencies.stairway;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.MissingRequiredFieldException;
import bio.terra.common.stairway.MonitoringHook;
import bio.terra.common.stairway.StairwayComponent;
import bio.terra.pipelines.common.utils.MdcHook;
import bio.terra.pipelines.stairway.CreateJobFlightMapKeys;
import bio.terra.stairway.Flight;
import bio.terra.stairway.FlightMap;
import java.util.UUID;
import javax.annotation.Nullable;
import org.apache.commons.lang3.StringUtils;

public class StairwayJobBuilder {
  private final StairwayJobService stairwayJobService;
  private final StairwayComponent stairwayComponent;
  private final MdcHook mdcHook;
  private final FlightMap jobParameterMap;
  private Class<? extends Flight> flightClass;
  private UUID jobId;
  @Nullable private String description;
  @Nullable private Object request;

  private String pipelineId;
  private String pipelineVersion;
  private String submittingUserId;
  private Object pipelineInputs;

  public StairwayJobBuilder(
      StairwayJobService stairwayJobService, StairwayComponent stairwayComponent, MdcHook mdcHook) {
    this.stairwayJobService = stairwayJobService;
    this.stairwayComponent = stairwayComponent;
    this.mdcHook = mdcHook;
    this.jobParameterMap = new FlightMap();
  }

  public StairwayJobBuilder flightClass(Class<? extends Flight> flightClass) {
    this.flightClass = flightClass;
    return this;
  }

  public StairwayJobBuilder jobId(UUID jobId) {
    this.jobId = jobId;
    return this;
  }

  public StairwayJobBuilder description(@Nullable String description) {
    this.description = description;
    return this;
  }

  public StairwayJobBuilder request(@Nullable Object request) {
    this.request = request;
    return this;
  }

  public StairwayJobBuilder pipelineId(@Nullable String pipelineId) {
    this.pipelineId = pipelineId;
    return this;
  }

  public StairwayJobBuilder pipelineVersion(@Nullable String pipelineVersion) {
    this.pipelineVersion = pipelineVersion;
    return this;
  }

  public StairwayJobBuilder submittingUserId(@Nullable String submittingUserId) {
    this.submittingUserId = submittingUserId;
    return this;
  }

  public StairwayJobBuilder pipelineInputs(@Nullable Object pipelineInputs) {
    this.pipelineInputs = pipelineInputs;
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
  public UUID submit() {
    populateInputParams();
    return stairwayJobService.submit(flightClass, jobParameterMap, jobId);
  }

  // Check the inputs, supply defaults and finalize the input parameter map
  private void populateInputParams() {
    if (flightClass == null) {
      throw new MissingRequiredFieldException("Missing flight class: flightClass");
    }

    // Always add the MDC logging and tracing span parameters for the mdc hook
    addParameter(MdcHook.MDC_FLIGHT_MAP_KEY, mdcHook.getSerializedCurrentContext());
    addParameter(
        MonitoringHook.SUBMISSION_SPAN_CONTEXT_MAP_KEY,
        MonitoringHook.serializeCurrentTracingContext());

    // Convert any other members that were set into parameters. However, if they were
    // explicitly added with addParameter during construction, we do not overwrite them.
    if (shouldInsert(StairwayJobMapKeys.DESCRIPTION, description)) {
      addParameter(StairwayJobMapKeys.DESCRIPTION.getKeyName(), description);
    }
    if (shouldInsert(StairwayJobMapKeys.REQUEST, request)) {
      addParameter(StairwayJobMapKeys.REQUEST.getKeyName(), request);
    }
    if (shouldInsert(CreateJobFlightMapKeys.PIPELINE_ID, pipelineId)) {
      addParameter(CreateJobFlightMapKeys.PIPELINE_ID, pipelineId);
    }
    if (shouldInsert(CreateJobFlightMapKeys.PIPELINE_VERSION, pipelineVersion)) {
      addParameter(CreateJobFlightMapKeys.PIPELINE_VERSION, pipelineVersion);
    }
    if (shouldInsert(CreateJobFlightMapKeys.SUBMITTING_USER_ID, submittingUserId)) {
      addParameter(CreateJobFlightMapKeys.SUBMITTING_USER_ID, submittingUserId);
    }
    if (shouldInsert(CreateJobFlightMapKeys.PIPELINE_INPUTS, pipelineInputs)) {
      addParameter(CreateJobFlightMapKeys.PIPELINE_INPUTS, pipelineInputs);
    }
  }

  private boolean shouldInsert(String mapKey, @Nullable Object value) {
    return (value != null && !jobParameterMap.containsKey(mapKey));
  }

  private boolean shouldInsert(StairwayJobMapKeys mapKey, @Nullable Object value) {
    return (value != null && !jobParameterMap.containsKey(mapKey.getKeyName()));
  }
}
