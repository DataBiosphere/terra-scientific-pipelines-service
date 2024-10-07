package bio.terra.pipelines.dependencies.stairway;

import java.util.Arrays;
import java.util.List;

public enum JobMapKeys {
  // parameters for all flight types
  DESCRIPTION("description"),
  USER_ID("user_id"),
  PIPELINE_NAME("pipeline_name"),
  STATUS_CODE("status_code"),
  RESPONSE("response"), // result or output of the job
  RESULT_PATH(
      "result_path"), // path to the result API endpoint for this job; only used for asynchronous

  // keys to determine which Stairway hooks to run
  DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK(
      "do_set_pipeline_run_status_failed_hook"), // whether to run the
  // StairwaySetPipelineRunStatusHook
  DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK(
      "do_increment_metrics_failed_counter_hook"); // whether to run the
  // StairwayFailedMetricsCounterHook

  private final String keyName;

  JobMapKeys(String keyName) {
    this.keyName = keyName;
  }

  public String getKeyName() {
    return keyName;
  }

  public static List<String> getRequiredKeys() {
    return Arrays.asList(JobMapKeys.USER_ID.getKeyName(), JobMapKeys.PIPELINE_NAME.getKeyName());
  }
}
