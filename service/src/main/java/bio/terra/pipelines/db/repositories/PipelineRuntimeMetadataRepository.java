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
   * Find the latest non-hidden runtime metadata entry for a given pipeline name
   *
   * @param pipelineKeyPrefix the pipeline name prefix (e.g., "array_imputation")
   * @return list of metadata entries matching the prefix
   */
  PipelineRuntimeMetadata findFirstByPipelineKeyStartingWithAndHiddenIsFalseOrderByPipelineKeyDesc(
      String pipelineKeyPrefix);

  /**
   * Find all runtime metadata entries with the given hidden status.
   *
   * @param hidden true to find hidden pipelines, false to find visible ones
   * @return list of matching metadata entries
   */
  List<PipelineRuntimeMetadata> findAllByHidden(boolean hidden);
}
