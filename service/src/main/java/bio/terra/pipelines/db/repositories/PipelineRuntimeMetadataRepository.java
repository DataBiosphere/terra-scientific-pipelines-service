package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.PipelineRuntimeMetadata;
import java.util.List;
import org.jetbrains.annotations.NotNull;
import org.springframework.data.repository.CrudRepository;

public interface PipelineRuntimeMetadataRepository
    extends CrudRepository<PipelineRuntimeMetadata, String> {

  @NotNull
  @Override
  List<PipelineRuntimeMetadata> findAll();
}
