package bio.terra.pipelines.dependencies.gcs;

import com.google.cloud.storage.Storage;
import com.google.cloud.storage.StorageOptions;
import org.springframework.stereotype.Component;

@Component
public class GcsClient {

  public GcsClient() {}

  public Storage getStorageService(String projectId) {
    // this will use application default credentials, which is what is used by
    // SamService.getTeaspoonsServiceAccountToken()
    return StorageOptions.newBuilder().setProjectId(projectId).build().getService();
  }
}
