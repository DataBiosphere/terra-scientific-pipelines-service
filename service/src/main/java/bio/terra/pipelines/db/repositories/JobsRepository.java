package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.Job;
import java.util.List;
import java.util.Optional;
import org.springframework.data.repository.CrudRepository;

public interface JobsRepository extends CrudRepository<Job, String> {
  List<Job> findAllByPipelineIdAndUserId(String pipelineId, String userId);

  Optional<Job> findJobByPipelineIdAndUserIdAndJobId(
      String pipelineId, String userId, String jobId);

  Optional<Job> findJobByJobId(String jobId);
}
