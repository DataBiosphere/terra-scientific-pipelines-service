package bio.terra.pipelines.dependencies.stairway;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.MissingRequiredFieldException;
import bio.terra.common.stairway.MonitoringHook;
import bio.terra.pipelines.common.utils.MdcHook;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.stairway.RunImputationJobFlightMapKeys;
import bio.terra.stairway.Flight;
import bio.terra.stairway.FlightMap;
import java.util.UUID;
import javax.annotation.Nullable;
import org.apache.commons.lang3.StringUtils;

public class JobBuilder {
  private final JobService jobService;
  private final MdcHook mdcHook;
  private final FlightMap jobParameterMap;
  private Class<? extends Flight> flightClass;
  private UUID jobId;
  @Nullable private String description;
  @Nullable private Object request;
  private String userId;

  private PipelinesEnum pipelineId;
  private String pipelineVersion;
  private Object pipelineInputs;

  public JobBuilder(JobService jobService, MdcHook mdcHook) {
    this.jobService = jobService;
    this.mdcHook = mdcHook;
    this.jobParameterMap = new FlightMap();
  }

  public JobBuilder flightClass(Class<? extends Flight> flightClass) {
    this.flightClass = flightClass;
    return this;
  }

  public JobBuilder jobId(UUID jobId) {
    this.jobId = jobId;
    return this;
  }

  public JobBuilder userId(String userId) {
    this.userId = userId;
    return this;
  }

  public JobBuilder pipelineId(PipelinesEnum pipelineId) {
    this.pipelineId = pipelineId;
    return this;
  }

  public JobBuilder description(@Nullable String description) {
    this.description = description;
    return this;
  }

  public JobBuilder request(@Nullable Object request) {
    this.request = request;
    return this;
  }

  public JobBuilder pipelineVersion(@Nullable String pipelineVersion) {
    this.pipelineVersion = pipelineVersion;
    return this;
  }

  public JobBuilder pipelineInputs(@Nullable Object pipelineInputs) {
    this.pipelineInputs = pipelineInputs;
    return this;
  }

  public JobBuilder addParameter(String keyName, @Nullable Object val) {
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
    return jobService.submit(flightClass, jobParameterMap, jobId);
  }

  // Check the inputs, supply defaults and finalize the input parameter map
  private void populateInputParams() {
    if (flightClass == null) {
      throw new MissingRequiredFieldException(
          "Missing required field for flight construction: flightClass");
    }

    if (jobId == null) {
      throw new MissingRequiredFieldException(
          "Missing required field for flight construction: jobId");
    }

    if (userId == null) {
      throw new MissingRequiredFieldException(
          "Missing required field for flight construction: userId");
    }

    if (pipelineId == null) {
      throw new MissingRequiredFieldException(
          "Missing required field for flight construction: pipelineId");
    }

    // Always add the MDC logging and tracing span parameters for the mdc hook
    addParameter(MdcHook.MDC_FLIGHT_MAP_KEY, mdcHook.getSerializedCurrentContext());
    addParameter(
        MonitoringHook.SUBMISSION_SPAN_CONTEXT_MAP_KEY,
        MonitoringHook.serializeCurrentTracingContext());

    // Convert any other members that were set into parameters. However, if they were
    // explicitly added with addParameter during construction, we do not overwrite them.
    if (shouldInsert(JobMapKeys.DESCRIPTION, description)) {
      addParameter(JobMapKeys.DESCRIPTION.getKeyName(), description);
    }
    if (shouldInsert(JobMapKeys.REQUEST, request)) {
      addParameter(JobMapKeys.REQUEST.getKeyName(), request);
    }
    if (shouldInsert(JobMapKeys.USER_ID, userId)) {
      addParameter(JobMapKeys.USER_ID.getKeyName(), userId);
    }
    if (shouldInsert(JobMapKeys.PIPELINE_ID, pipelineId)) {
      addParameter(JobMapKeys.PIPELINE_ID.getKeyName(), pipelineId);
    }
    if (shouldInsert(RunImputationJobFlightMapKeys.PIPELINE_VERSION, pipelineVersion)) {
      addParameter(RunImputationJobFlightMapKeys.PIPELINE_VERSION, pipelineVersion);
    }
    if (shouldInsert(RunImputationJobFlightMapKeys.PIPELINE_INPUTS, pipelineInputs)) {
      addParameter(RunImputationJobFlightMapKeys.PIPELINE_INPUTS, pipelineInputs);
    }
  }

  private boolean shouldInsert(String mapKey, @Nullable Object value) {
    return (value != null && !jobParameterMap.containsKey(mapKey));
  }

  private boolean shouldInsert(JobMapKeys mapKey, @Nullable Object value) {
    return (value != null && !jobParameterMap.containsKey(mapKey.getKeyName()));
  }
}
