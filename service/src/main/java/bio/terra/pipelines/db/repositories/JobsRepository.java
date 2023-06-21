package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.DbJob;
import java.util.List;
import java.util.Optional;
import org.springframework.data.repository.CrudRepository;

public interface JobsRepository extends CrudRepository<DbJob, Long> {
  List<DbJob> findAllByPipelineIdAndUserId(String pipelineId, String userId);

  Optional<DbJob> findJobByPipelineIdAndUserIdAndJobId(
      String pipelineId, String userId, String jobId);

  Optional<DbJob> findJobByJobId(String jobId);
}
