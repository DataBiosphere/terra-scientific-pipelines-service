package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import java.util.List;
import org.jetbrains.annotations.NotNull;
import org.springframework.data.repository.CrudRepository;

public interface PipelinesRepository extends CrudRepository<Pipeline, Long> {

  @NotNull
  @Override
  List<Pipeline> findAll();

  Boolean existsByName(PipelinesEnum name);

  Pipeline findByName(PipelinesEnum name);
}
