package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.DataDelivery;
import java.util.Optional;
import org.springframework.data.repository.CrudRepository;

public interface DataDeliveryRepository extends CrudRepository<DataDelivery, Long> {
  Optional<DataDelivery> findFirstByPipelineRunIdOrderByCreatedDesc(Long pipelineRunId);
}
