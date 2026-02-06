package bio.terra.pipelines.dependencies.gcs;

import com.google.auth.oauth2.AccessToken;
import com.google.auth.oauth2.GoogleCredentials;
import com.google.cloud.storage.Storage;
import com.google.cloud.storage.StorageOptions;
import org.springframework.stereotype.Component;

@Component
public class GcsClient {
  public Storage getStorageService(String projectId) {
    // this will use application default credentials, which is what is used by
    // SamService.getTeaspoonsServiceAccountToken()
    return StorageOptions.newBuilder().setProjectId(projectId).build().getService();
  }

  public Storage getStorageService(String projectId, String bearerToken) {
    if (bearerToken == null || bearerToken.isEmpty()) {
      return getStorageService(projectId);
    }
    GoogleCredentials userCredentials = userCredentialsFromBearerToken(bearerToken);
    return StorageOptions.newBuilder()
        .setProjectId(projectId)
        .setCredentials(userCredentials)
        .build()
        .getService();
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
