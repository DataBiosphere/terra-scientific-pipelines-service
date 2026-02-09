package bio.terra.pipelines.dependencies.gcs;

import com.google.auth.oauth2.AccessToken;
import com.google.auth.oauth2.GoogleCredentials;
import com.google.cloud.storage.Storage;
import com.google.cloud.storage.StorageOptions;
import org.springframework.stereotype.Component;

@Component
public class GcsClient {
  /**
   * Get a Storage service instance using application default credentials and the specified project
   * ID.
   */
  public Storage getStorageServiceWithProject(String projectId) {
    // this will use application default credentials, which is what is used by
    // SamService.getTeaspoonsServiceAccountToken()
    return StorageOptions.newBuilder().setProjectId(projectId).build().getService();
  }

  /**
   * Get a Storage service instance using the provided bearer token. If the token is null or empty,
   * use application default credentials.
   *
   * @param bearerToken string value of a token, or null/empty to use application default
   *     credentials
   */
  public Storage getStorageService(String bearerToken) {
    if (bearerToken == null || bearerToken.isEmpty()) {
      return StorageOptions.newBuilder().build().getService();
    }
    GoogleCredentials userCredentials = userCredentialsFromBearerToken(bearerToken);
    return StorageOptions.newBuilder().setCredentials(userCredentials).build().getService();
  }

  /**
   * Create GoogleCredentials from a bearer token. The token is expected to be a user access token.
   *
   * @param bearerToken string value of the token
   * @return GoogleCredentials object
   */
  private GoogleCredentials userCredentialsFromBearerToken(String bearerToken) {
    AccessToken accessToken = new AccessToken(bearerToken, null);
    return GoogleCredentials.create(accessToken);
  }
}
