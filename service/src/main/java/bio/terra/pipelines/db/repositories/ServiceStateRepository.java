package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.ServiceState;
import org.springframework.data.repository.CrudRepository;

public interface ServiceStateRepository extends CrudRepository<ServiceState, Integer> {}
