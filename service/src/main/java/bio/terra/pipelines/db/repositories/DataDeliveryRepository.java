package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.DataDelivery;
import java.util.List;
import java.util.Optional;
import java.util.UUID;
import org.springframework.data.repository.CrudRepository;

public interface DataDeliveryRepository extends CrudRepository<DataDelivery, Long> {
  Optional<DataDelivery> findByJobId(UUID jobId);

  List<DataDelivery> findAllByPipelineRunId(Long pipelineRunId);

  Optional<DataDelivery> findFirstByPipelineRunIdOrderByCreatedDesc(Long pipelineRunId);

  List<DataDelivery> findAllByStatus(String status);

  boolean existsByJobId(UUID jobId);
}
