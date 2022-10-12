package bio.terra.pipelines.db;

import bio.terra.pipelines.common.model.PolicyInput;

/** Record to hold one attribute of an attribute set and its set id */
public record DbAttribute(String setId, PolicyInput policyInput) {}
