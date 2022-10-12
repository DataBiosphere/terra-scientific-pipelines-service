package bio.terra.pipelines.service.pipeline;

import bio.terra.pipelines.common.exception.PolicyNotImplementedException;
import bio.terra.pipelines.common.model.PolicyInputs;
import bio.terra.pipelines.service.pipeline.model.CloudPlatform;
import bio.terra.pipelines.service.pipeline.model.RegionConstraintResult;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class PolicyService {

  @Autowired
  public PolicyService() {}

  public RegionConstraintResult evaluateRegionConstraint(
      PolicyInputs policyInputs, CloudPlatform cloudPlatform, String regionRequest) {
    throw new PolicyNotImplementedException("Not implemented");
  }
}
