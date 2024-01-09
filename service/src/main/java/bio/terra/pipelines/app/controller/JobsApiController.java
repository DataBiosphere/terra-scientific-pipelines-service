package bio.terra.pipelines.app.controller;

import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.dependencies.stairway.model.EnumeratedJobs;
import bio.terra.pipelines.generated.api.JobsApi;
import bio.terra.pipelines.generated.model.*;
import bio.terra.stairway.FlightState;
import io.swagger.annotations.Api;
import java.util.UUID;
import javax.servlet.http.HttpServletRequest;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.PathVariable;

/** Jobs controller */
@Controller
@Api(tags = {"jobs"})
public class JobsApiController implements JobsApi {
  private final SamConfiguration samConfiguration;
  private final SamUserFactory samUserFactory;
  private final HttpServletRequest request;
  private final StairwayJobService stairwayJobService;

  @Autowired
  public JobsApiController(
      SamConfiguration samConfiguration,
      SamUserFactory samUserFactory,
      HttpServletRequest request,
      StairwayJobService stairwayJobService) {
    this.samConfiguration = samConfiguration;
    this.samUserFactory = samUserFactory;
    this.request = request;
    this.stairwayJobService = stairwayJobService;
  }

  private static final Logger logger = LoggerFactory.getLogger(JobsApiController.class);

  private SamUser getAuthenticatedInfo() {
    return samUserFactory.from(request, samConfiguration.baseUri());
  }

  // -- Jobs --

  @Override
  public ResponseEntity<ApiJobReport> getJob(@PathVariable("jobId") UUID jobId) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    FlightState flightState = stairwayJobService.retrieveJob(jobId, userId);
    ApiJobReport result = JobApiUtils.mapFlightStateToApiJobReport(flightState);
    return new ResponseEntity<>(result, HttpStatus.OK);
  }

  @Override
  public ResponseEntity<ApiGetJobsResponse> getAllJobs(Integer limit, String pageToken) {
    final SamUser userRequest = getAuthenticatedInfo();
    String userId = userRequest.getSubjectId();
    EnumeratedJobs enumeratedJobs =
        stairwayJobService.enumerateJobs(userId, limit, pageToken, null);
    ApiGetJobsResponse result = JobApiUtils.mapEnumeratedJobsToApi(enumeratedJobs);
    return new ResponseEntity<>(result, HttpStatus.OK);
  }
}
