package bio.terra.pipelines.stairway.flights.imputation;

public class ImputationJobMapKeys {
  public static final String PIPELINE_INPUT_DEFINITIONS = "pipeline_input_definitions";
  public static final String PIPELINE_OUTPUT_DEFINITIONS = "pipeline_output_definitions";
  public static final String USER_PROVIDED_PIPELINE_INPUTS = "user_provided_pipeline_inputs";
  public static final String ALL_PIPELINE_INPUTS = "all_pipeline_inputs";
  public static final String CONTROL_WORKSPACE_STORAGE_CONTAINER_NAME =
      "control_workspace_storage_container_name";
  public static final String CONTROL_WORKSPACE_STORAGE_CONTAINER_PROTOCOL =
      "control_workspace_storage_container_protocol";
  public static final String WDL_METHOD_NAME = "wdl_method_name";
  public static final String WDL_METHOD_VERSION = "wdl_method_version";
  public static final String PIPELINE_RUN_OUTPUTS = "pipeline_run_outputs";

  // GCP specific keys
  public static final String CONTROL_WORKSPACE_BILLING_PROJECT =
      "control_workspace_billing_project";
  public static final String CONTROL_WORKSPACE_NAME = "control_workspace_name";
  public static final String SUBMISSION_ID = "submission_id";
  public static final String QUOTA_SUBMISSION_ID = "quota_submission_id";
  public static final String RAW_QUOTA_CONSUMED = "raw_quota_consumed";
  public static final String EFFECTIVE_QUOTA_CONSUMED = "effective_quota_consumed";

  // Azure specific keys
  public static final String CONTROL_WORKSPACE_ID = "control_workspace_id";
  public static final String CBAS_URI = "cbas_uri";
  public static final String WDS_URI = "wds_uri";
  public static final String RUN_SET_ID = "run_set_id";

  ImputationJobMapKeys() {
    throw new IllegalStateException("Attempted to instantiate utility class ImputationJobMapKeys");
  }
}
