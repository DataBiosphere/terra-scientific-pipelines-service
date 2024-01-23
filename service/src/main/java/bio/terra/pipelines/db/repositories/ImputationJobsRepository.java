package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.ImputationJob;
import java.util.List;
import java.util.Optional;
import java.util.UUID;
import org.springframework.data.repository.CrudRepository;

public interface ImputationJobsRepository extends CrudRepository<ImputationJob, Long> {
  List<ImputationJob> findAllByUserId(String userId);

  Optional<ImputationJob> findJobByJobIdAndUserId(UUID jobId, String userId);
}
