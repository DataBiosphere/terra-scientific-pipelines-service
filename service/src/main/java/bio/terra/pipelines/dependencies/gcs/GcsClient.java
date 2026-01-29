package bio.terra.pipelines.dependencies.gcs;

import bio.terra.common.iam.BearerToken;
import com.google.auth.oauth2.AccessToken;
import com.google.auth.oauth2.GoogleCredentials;
import com.google.cloud.storage.Storage;
import com.google.cloud.storage.StorageOptions;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

@Component
public class GcsClient {
  private static final Logger logger = LoggerFactory.getLogger(GcsClient.class);

  public Storage getStorageService(String projectId) {
    // this will use application default credentials, which is what is used by
    // SamService.getTeaspoonsServiceAccountToken()
    return StorageOptions.newBuilder().setProjectId(projectId).build().getService();
  }

  public Storage getUserStorageService(String projectId, BearerToken bearerToken) {
    //    logger.info("using bearer token {}", bearerToken.getToken());
    GoogleCredentials userCredentials = userCredentialsFromBearerToken(bearerToken);
    //    logger.info("created user credentials with access token {}",
    // userCredentials.getAccessToken());
    return StorageOptions.newBuilder()
        .setProjectId(projectId)
        .setCredentials(userCredentials)
        .build()
        .getService();
  }

  private GoogleCredentials userCredentialsFromBearerToken(BearerToken bearerToken) {
    //    Date expirationTime = new Date(System.currentTimeMillis() + 3600 * 1000); // expires in 1
    // hour

    // 1. Create an AccessToken object
    AccessToken accessToken = new AccessToken(bearerToken.getToken(), null);

    // 2. Create GoogleCredentials from the AccessToken
    return GoogleCredentials.create(accessToken);
  }
}
