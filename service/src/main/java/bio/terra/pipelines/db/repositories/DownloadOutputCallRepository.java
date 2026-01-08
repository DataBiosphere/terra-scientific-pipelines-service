package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.DownloadOutputCall;
import java.util.Optional;
import java.util.UUID;
import org.springframework.data.repository.CrudRepository;

public interface DownloadOutputCallRepository extends CrudRepository<DownloadOutputCall, Long> {
  Optional<DownloadOutputCall> findByJobId(UUID jobId);
}
