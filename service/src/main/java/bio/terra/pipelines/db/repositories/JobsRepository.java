package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.Job;
import java.util.List;
import java.util.Optional;
import java.util.UUID;
import org.springframework.data.repository.CrudRepository;

public interface JobsRepository extends CrudRepository<Job, Long> {
  List<Job> findAllByUserId(String userId);

  Optional<Job> findJobByJobIdAndUserId(UUID jobId, String userId);
}
