package bio.terra.pipelines.dependencies.gcs;

import com.google.auth.oauth2.AccessToken;
import com.google.auth.oauth2.GoogleCredentials;
import com.google.cloud.storage.Storage;
import com.google.cloud.storage.StorageOptions;
import org.springframework.stereotype.Component;

@Component
public class GcsClient {
  public Storage getStorageServiceWithProject(String projectId) {
    // this will use application default credentials, which is what is used by
    // SamService.getTeaspoonsServiceAccountToken()
    return StorageOptions.newBuilder().setProjectId(projectId).build().getService();
  }

  /**
   * Get a Storage service instance using the provided bearer token. If the token is null or empty,
   * use application default credentials.
   */
  public Storage getStorageService(String bearerToken) {
    if (bearerToken == null || bearerToken.isEmpty()) {
      // this will use application default credentials
      return StorageOptions.newBuilder().build().getService();
    }
    GoogleCredentials userCredentials = userCredentialsFromBearerToken(bearerToken);
    return StorageOptions.newBuilder().setCredentials(userCredentials).build().getService();
  }

  private GoogleCredentials userCredentialsFromBearerToken(String bearerToken) {
    //    Date expirationTime = new Date(System.currentTimeMillis() + 3600 * 1000); // expires in 1
    // hour

    // 1. Create an AccessToken object
    AccessToken accessToken = new AccessToken(bearerToken, null);

    // 2. Create GoogleCredentials from the AccessToken
    return GoogleCredentials.create(accessToken);
  }
}
