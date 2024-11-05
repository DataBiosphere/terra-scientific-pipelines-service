package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.UserQuota;
import java.util.ArrayList;
import java.util.Optional;
import org.springframework.data.repository.CrudRepository;

public interface UserQuotasRepository extends CrudRepository<UserQuota, Long> {
  Optional<UserQuota> findByUserIdAndPipelineName(String userId, PipelinesEnum pipelineName);

  ArrayList<UserQuota> findByUserId(String userId);
}
