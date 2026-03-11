package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.DataDelivery;
import java.util.List;
import java.util.Optional;
import org.springframework.data.repository.CrudRepository;

public interface DataDeliveryRepository extends CrudRepository<DataDelivery, Long> {
  List<DataDelivery> findAllByPipelineRunId(Long pipelineRunId);

  Optional<DataDelivery> findFirstByPipelineRunIdOrderByCreatedDesc(Long pipelineRunId);
}
