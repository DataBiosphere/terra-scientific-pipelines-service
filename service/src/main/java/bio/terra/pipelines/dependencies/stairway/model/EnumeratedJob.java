package bio.terra.pipelines.dependencies.stairway.model;

import bio.terra.stairway.FlightState;

/** Class to store a Stairway job result that translates nicely into API responses */
public class EnumeratedJob {
  private FlightState flightState;
  private String jobDescription;

  public FlightState getFlightState() {
    return flightState;
  }

  public EnumeratedJob flightState(FlightState flightState) {
    this.flightState = flightState;
    return this;
  }

  public String getJobDescription() {
    return jobDescription;
  }

  public EnumeratedJob jobDescription(String jobDescription) {
    this.jobDescription = jobDescription;
    return this;
  }
}
