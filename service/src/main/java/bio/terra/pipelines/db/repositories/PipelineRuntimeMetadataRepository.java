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

  /**
   * Find all runtime metadata entries for a given pipeline name.
   *
   * @param pipelineKeyPrefix the pipeline name prefix (e.g., "array_imputation_")
   * @return list of metadata entries matching the prefix
   */
  List<PipelineRuntimeMetadata> findAllByPipelineKeyStartingWith(String pipelineKeyPrefix);
}
