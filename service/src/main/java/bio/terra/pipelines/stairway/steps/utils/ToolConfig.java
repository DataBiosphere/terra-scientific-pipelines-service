package bio.terra.pipelines.stairway.steps.utils;

import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.db.entities.PipelineOutputDefinition;
import java.util.List;

public record ToolConfig(
    String methodName,
    String methodVersion,
    List<PipelineInputDefinition> inputDefinitions,
    List<PipelineOutputDefinition> outputDefinitions,
    boolean callCache,
    boolean deleteIntermediateOutputFiles,
    boolean useReferenceDisks,
    Long pollingIntervalSeconds) {
  // class to store relevant information for submitting a method/workflow

}
