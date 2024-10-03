package bio.terra.pipelines.stairway.imputation;

import bio.terra.stairway.Flight;
import bio.terra.stairway.FlightMap;

/**
 * Class to be implemented by any flight class that deals with PipelineRuns and therefore should
 * have relevant hooks (StairwaySetPipelineRunStatusHook and StairwayFailedMetricsCounterHook)
 * applied.
 */
public class PipelineRunTypeFlight extends Flight {
  public PipelineRunTypeFlight(FlightMap inputParameters, Object applicationContext) {
    super(inputParameters, applicationContext);
  }
}
