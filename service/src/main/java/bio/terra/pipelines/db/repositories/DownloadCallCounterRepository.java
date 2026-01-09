package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.DownloadCallCount;
import java.util.Optional;
import java.util.UUID;
import org.springframework.data.repository.CrudRepository;

public interface DownloadCallCounterRepository extends CrudRepository<DownloadCallCount, Long> {
  Optional<DownloadCallCount> findByJobId(UUID jobId);
}
