package bio.terra.pipelines.service;

import bio.terra.pipelines.app.configuration.internal.PipelineCatalogConfigurations;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineQuota;
import bio.terra.pipelines.db.entities.PipelineRuntimeMetadata;
import bio.terra.pipelines.model.Pipeline;
import java.util.List;
import java.util.Optional;

/** Domain-facing catalog contract for pipeline definitions and assembly. */
public interface PipelineCatalogRepository {
  Optional<PipelineCatalogConfigurations.PipelineDefinition> getDefinition(
      PipelinesEnum pipelineName, Integer version);

  List<PipelineCatalogService.ConfiguredPipelineVersion> getConfiguredPipelineVersions();

  Pipeline hydratePipeline(PipelineRuntimeMetadata dbPipeline);

  Pipeline buildPipeline(
      PipelinesEnum pipelineName, Integer version, PipelineRuntimeMetadata dbPipelineOrNull);

  PipelineQuota getPipelineQuota(PipelinesEnum pipelineName);
}
