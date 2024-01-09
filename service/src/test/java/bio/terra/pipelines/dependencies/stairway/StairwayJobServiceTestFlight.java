package bio.terra.pipelines.dependencies.stairway;

import bio.terra.stairway.Flight;
import bio.terra.stairway.FlightMap;

public class StairwayJobServiceTestFlight extends Flight {

  public StairwayJobServiceTestFlight(FlightMap inputParameters, Object applicationContext) {
    super(inputParameters, applicationContext);

    // Just one step for this test
    addStep(new StairwayJobServiceTestStep());
  }
}
