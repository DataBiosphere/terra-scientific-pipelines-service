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
   * @param pipelineName the pipeline name (e.g., "array_imputation")
   * @return list of metadata entries matching the name
   */
  PipelineRuntimeMetadata findFirstByPipelineNameAndHiddenIsFalseOrderByPipelineKeyDesc(
      String pipelineName);
}
