package bio.terra.pipelines.stairway.flights.datadelivery;

public class DataDeliveryJobMapKeys {
  public static final String DESTINATION_GCS_PATH = "destination_gcs_path";
  public static final String PIPELINE_RUN_ID = "pipeline_run_id";

  DataDeliveryJobMapKeys() {
    throw new IllegalStateException(
        "Attempted to instantiate utility class DataDeliveryJobMapKeys");
  }
}
