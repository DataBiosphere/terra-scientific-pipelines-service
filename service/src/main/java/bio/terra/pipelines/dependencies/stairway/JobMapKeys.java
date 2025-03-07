package bio.terra.pipelines.dependencies.stairway;

import java.util.Arrays;
import java.util.List;

public class JobMapKeys {
  // parameters for all flight types
  public static final String DESCRIPTION = "description";
  public static final String USER_ID = "user_id";
  public static final String PIPELINE_NAME = "pipeline_name";
  public static final String PIPELINE_ID = "pipeline_id";
  public static final String STATUS_CODE = "status_code";
  public static final String RESPONSE = "response"; // result or output of the job
  // domain name for the service, used to generate the URL for the result api endpoint
  public static final String DOMAIN_NAME = "domain_name";

  // keys to determine which Stairway hooks to run
  public static final String DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK =
      "do_set_pipeline_run_status_failed_hook";
  public static final String DO_INCREMENT_METRICS_FAILED_COUNTER_HOOK =
      "do_increment_metrics_failed_counter_hook";
  public static final String DO_SEND_JOB_FAILURE_NOTIFICATION_HOOK =
      "do_send_job_failure_notification_hook";

  JobMapKeys() {
    throw new IllegalStateException("Attempted to instantiate utility class JobMapKeys");
  }

  public static List<String> getRequiredKeys() {
    return Arrays.asList(USER_ID, PIPELINE_NAME, PIPELINE_ID);
  }
}
