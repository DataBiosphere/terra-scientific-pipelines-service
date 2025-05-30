package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.PipelineRun;
import java.util.List;
import java.util.Optional;
import java.util.UUID;
import org.springframework.data.jpa.repository.JpaSpecificationExecutor;
import org.springframework.data.repository.CrudRepository;

public interface PipelineRunsRepository
    extends CrudRepository<PipelineRun, Long>, JpaSpecificationExecutor<PipelineRun> {
  List<PipelineRun> findAllByUserId(String userId);

  boolean existsByJobId(UUID jobId);

  Optional<PipelineRun> findByJobIdAndUserId(UUID jobId, String userId);

  boolean existsByIdGreaterThan(Long id);

  int countByUserId(String userId);
}
