package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.Job;
import java.util.List;
import java.util.Optional;
import java.util.UUID;
import org.springframework.data.repository.CrudRepository;

public interface JobsRepository extends CrudRepository<Job, Long> {
  List<Job> findAllByPipelineIdAndUserId(String pipelineId, String userId);

  Optional<Job> findJobByPipelineIdAndUserIdAndJobId(String pipelineId, String userId, UUID jobId);
}
