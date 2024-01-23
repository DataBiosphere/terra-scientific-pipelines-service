package bio.terra.pipelines.dependencies.stairway;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.MissingRequiredFieldException;
import bio.terra.common.stairway.MonitoringHook;
import bio.terra.pipelines.common.utils.MdcHook;
import bio.terra.stairway.Flight;
import bio.terra.stairway.FlightMap;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.UUID;
import javax.annotation.Nullable;
import org.apache.commons.lang3.StringUtils;

public class JobBuilder {
  private final JobService jobService;
  private final MdcHook mdcHook;
  private final FlightMap jobParameterMap;
  private Class<? extends Flight> flightClass;
  private UUID jobId;

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
    validateRequiredInputs();
    return jobService.submit(flightClass, jobParameterMap, jobId);
  }

  // Supply defaults and finalize the input parameter map
  private void populateInputParams() {
    // Always add the MDC logging and tracing span parameters for the mdc hook
    addParameter(MdcHook.MDC_FLIGHT_MAP_KEY, mdcHook.getSerializedCurrentContext());
    addParameter(
        MonitoringHook.SUBMISSION_SPAN_CONTEXT_MAP_KEY,
        MonitoringHook.serializeCurrentTracingContext());
  }

  private void validateRequiredInputs() {
    if (flightClass == null) {
      throw new MissingRequiredFieldException(
          "Missing required field for flight construction: flightClass");
    }

    if (jobId == null) {
      throw new MissingRequiredFieldException(
          "Missing required field for flight construction: jobId");
    }

    List<String> missingFields = new ArrayList<>();
    for (String keyName : JobMapKeys.getRequiredKeys()) {
      if (!jobParameterMap.containsKey(keyName)
          || Objects.equals(
              jobParameterMap.getRaw(keyName), "null") // getRaw stringifies the result
          || Objects.requireNonNull(jobParameterMap.getRaw(keyName)).equals("\"\"")) {
        missingFields.add(keyName);
      }
      if (!missingFields.isEmpty()) {
        throw new MissingRequiredFieldException(
            "Missing required field(s) for flight construction: "
                + String.join(", ", missingFields));
      }
    }
  }
}
