package bio.terra.pipelines.dependencies.sam;

import bio.terra.common.iam.BearerToken;
import bio.terra.common.sam.SamRetry;
import bio.terra.common.sam.exception.SamExceptionFactory;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.generated.model.ApiSystemStatusSystems;
import org.broadinstitute.dsde.workbench.client.sam.ApiException;
import org.broadinstitute.dsde.workbench.client.sam.model.SystemStatus;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class SamService implements HealthCheck {
  private static final Logger logger = LoggerFactory.getLogger(SamService.class);
  private final SamClient samClient;

  @Autowired
  public SamService(SamClient samClient) {
    this.samClient = samClient;
  }

  public boolean getAction(
      String resourceType, String resourceId, String action, BearerToken bearerToken) {
    try {
      return SamRetry.retry(
          () ->
              samClient
                  .resourcesApi(bearerToken.getToken())
                  .resourcePermissionV2(resourceType, resourceId, action));
    } catch (ApiException e) {
      throw SamExceptionFactory.create(e);
    } catch (InterruptedException e) {
      Thread.currentThread().interrupt();
      throw SamExceptionFactory.create("Sam retry interrupted", e);
    }
  }

  public Result checkHealth() {
    // No access token needed since this is an unauthenticated API.
    try {
      SystemStatus result = samClient.statusApi().getSystemStatus();
      return new Result(result.getOk(), result.toString());
    } catch (ApiException e) {
      String errorMsg = "Sam status check failed";
      logger.error(errorMsg, e);
      return new Result(false, e.getMessage());
    }
  }

  public ApiSystemStatusSystems checkHealthApiSystemStatus() {
    Result healthResult = checkHealth();
    return new ApiSystemStatusSystems()
        .ok(healthResult.isOk())
        .addMessagesItem(healthResult.message());
  }
}
