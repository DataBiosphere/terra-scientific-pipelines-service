package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineRuntimeMetadata;
import java.util.List;
import org.jetbrains.annotations.NotNull;
import org.springframework.data.repository.CrudRepository;

public interface PipelineRuntimeMetadataRepository
    extends CrudRepository<PipelineRuntimeMetadata, Long> {

  @NotNull
  @Override
  List<PipelineRuntimeMetadata> findAll();

  List<PipelineRuntimeMetadata> findAllByOrderByNameAscVersionDesc();

  List<PipelineRuntimeMetadata> findAllByHiddenIsFalseOrderByNameAscVersionDesc();

  Boolean existsByNameAndHiddenIsFalse(PipelinesEnum name);

  PipelineRuntimeMetadata findByNameAndHiddenIsFalse(PipelinesEnum name);

  PipelineRuntimeMetadata findByNameAndVersion(PipelinesEnum name, Integer pipelineVersion);

  PipelineRuntimeMetadata findByNameAndVersionAndHiddenIsFalse(
      PipelinesEnum name, Integer pipelineVersion);

  PipelineRuntimeMetadata findFirstByNameAndHiddenIsFalseOrderByVersionDesc(PipelinesEnum name);
}
