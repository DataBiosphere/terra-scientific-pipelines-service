package bio.terra.pipelines.testutils;

import bio.terra.pipelines.dependencies.stairway.JobServiceTestFlight;
import bio.terra.pipelines.dependencies.stairway.JobServiceTestStep;
import bio.terra.stairway.*;
import java.util.List;

/**
 * A flight context implementation for use in unit tests to avoid running Stairway flights, with
 * additional helper methods for modification.
 */
public class TestFlightContext implements FlightContext {

  private String flightId = TestUtils.TEST_NEW_UUID.toString();
  private String flightClassName = JobServiceTestFlight.class.getName();
  private FlightMap inputParameters = new FlightMap();
  private FlightMap workingMap = new FlightMap();
  private int stepIndex = 0;
  private FlightStatus flightStatus = FlightStatus.QUEUED;
  private Direction direction = Direction.DO;
  private String stepClassName = JobServiceTestStep.class.getName();
  private StepResult result = new StepResult(StepStatus.STEP_RESULT_SUCCESS);

  @Override
  public Object getApplicationContext() {
    return null;
  }

  @Override
  public String getFlightId() {
    return flightId;
  }

  public TestFlightContext flightId(String flightId) {
    this.flightId = flightId;
    return this;
  }

  @Override
  public String getFlightClassName() {
    return flightClassName;
  }

  @Override
  public FlightMap getInputParameters() {
    return inputParameters;
  }

  public TestFlightContext inputParameters(FlightMap inputParameters) {
    this.inputParameters = inputParameters;
    return this;
  }

  @Override
  public FlightMap getWorkingMap() {
    return workingMap;
  }

  public TestFlightContext workingMap(FlightMap workingMap) {
    this.workingMap = workingMap;
    return this;
  }

  @Override
  public int getStepIndex() {
    return stepIndex;
  }

  @Override
  public FlightStatus getFlightStatus() {
    return flightStatus;
  }

  public TestFlightContext flightStatus(FlightStatus flightStatus) {
    this.flightStatus = flightStatus;
    return this;
  }

  @Override
  public boolean isRerun() {
    return false;
  }

  @Override
  public Direction getDirection() {
    return direction;
  }

  @Override
  public StepResult getResult() {
    return result;
  }

  public TestFlightContext result(StepResult result) {
    this.result = result;
    return this;
  }

  @Override
  public Stairway getStairway() {
    return null;
  }

  @Override
  public List<String> getStepClassNames() {
    return null;
  }

  @Override
  public String getStepClassName() {
    return stepClassName;
  }

  @Override
  public String prettyStepState() {
    return null;
  }

  @Override
  public String flightDesc() {
    return null;
  }

  @Override
  public ProgressMeter getProgressMeter(String name) {
    return null;
  }

  @Override
  public void setProgressMeter(String name, long v1, long v2) throws InterruptedException {
    // no-op
  }
}
