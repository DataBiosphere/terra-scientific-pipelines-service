package bio.terra.pipelines.dependencies.sam;

import bio.terra.common.exception.ForbiddenException;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.common.iam.BearerToken;
import bio.terra.common.iam.SamUser;
import bio.terra.common.sam.SamRetry;
import bio.terra.common.sam.exception.SamExceptionFactory;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.generated.model.ApiSystemStatusSystems;
import com.google.auth.oauth2.GoogleCredentials;
import java.io.IOException;
import java.util.Set;
import org.broadinstitute.dsde.workbench.client.sam.ApiException;
import org.broadinstitute.dsde.workbench.client.sam.model.SystemStatus;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

/** Encapsulates logic for interacting with SAM */
@Component
public class SamService implements HealthCheck {
  private static final Set<String> SAM_OAUTH_SCOPES = Set.of("openid", "email", "profile");
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

  public String getTeaspoonsServiceAccountToken() {
    try {
      GoogleCredentials creds =
          GoogleCredentials.getApplicationDefault().createScoped(SAM_OAUTH_SCOPES);
      creds.refreshIfExpired();
      return creds.getAccessToken().getTokenValue();
    } catch (IOException e) {
      throw new InternalServerErrorException(
          "Internal server error retrieving Teaspoons credentials", e);
    }
  }

  /**
   * Wrapper around isAdmin which throws an appropriate exception if a user does not have admin
   * access.
   *
   * @param authenticatedUser Authenticated Sam user whose permissions are being checked
   */
  public void checkAdminAuthz(SamUser authenticatedUser) {
    boolean isAuthorized =
        isAdmin(authenticatedUser.getEmail(), authenticatedUser.getBearerToken().getToken());
    if (!isAuthorized)
      throw new ForbiddenException(
          String.format(
              "User %s is not authorized to perform admin action", authenticatedUser.getEmail()));
    else logger.info("User {} is an authorized admin", authenticatedUser.getEmail());
  }

  private boolean isAdmin(String userEmail, String accessToken) {
    try {
      // If the user can successfully call sam admin api, the user has terra level admin access.
      SamRetry.retry(() -> samClient.adminApi(accessToken).adminGetUserByEmail(userEmail));
      return true;
    } catch (ApiException apiException) {
      logger.info(
          "Error checking admin permission in Sam. This is expected if requester is not SAM admin.",
          apiException);
      return false;
    } catch (InterruptedException e) {
      Thread.currentThread().interrupt();
      throw SamExceptionFactory.create("Sam retry interrupted", e);
    }
  }
}
