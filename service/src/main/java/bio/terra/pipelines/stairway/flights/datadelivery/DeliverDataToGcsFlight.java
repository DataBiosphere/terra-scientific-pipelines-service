package bio.terra.pipelines.stairway.flights.datadelivery;

import bio.terra.stairway.Flight;
import bio.terra.stairway.FlightMap;

public class DeliverDataToGcsFlight extends Flight {
  public DeliverDataToGcsFlight(FlightMap inputParameters, Object beanBag) {
    super(inputParameters, beanBag);
  }
}
