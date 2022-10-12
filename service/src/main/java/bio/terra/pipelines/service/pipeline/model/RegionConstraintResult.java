package bio.terra.pipelines.service.pipeline.model;

import java.util.List;

public record RegionConstraintResult(
    String error, String message, String region, List<String> sites) {}
