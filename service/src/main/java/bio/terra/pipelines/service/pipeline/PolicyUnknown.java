package bio.terra.pipelines.service.pipeline;

import bio.terra.pipelines.common.model.PolicyInput;
import bio.terra.pipelines.common.model.PolicyName;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

/** When we process a policy that is not known to TPS, then we use this combiner. */
public class PolicyUnknown implements PolicyBase {
  private final PolicyName policyName;

  public PolicyUnknown(PolicyName policyName) {
    this.policyName = policyName;
  }

  @Override
  public PolicyName getPolicyName() {
    return policyName;
  }

  /**
   * Unknown combining uses a simplistic algorithm. If the policy inputs have no additional data,
   * then we return a policy input with the name and empty additional data. That has the effect of
   * treating the policy as a flag; e.g., if this has PHI and that has PHI, then the result has PHI.
   *
   * <p>It the policy inputs have additional data, the additional data must match exactly.
   * Otherwise, it is treated as a conflict.
   *
   * @param dependent policy input
   * @param source policy input
   * @return policy input or null, if there is a conflict
   */
  @Override
  public PolicyInput combine(PolicyInput dependent, PolicyInput source) {
    Multimap<String, String> dependentData = dependent.getAdditionalData();
    Multimap<String, String> sourceData = source.getAdditionalData();
    Multimap<String, String> newData = ArrayListMultimap.create();

    // "Flag" mode - no possible conflict; empty data
    if (dependentData.size() == 0 && sourceData.size() == 0) {
      return new PolicyInput(dependent.getPolicyName(), newData);
    }

    // "Exact match" test - conflict if data is not equal
    // We make a copy of the data, since (theoretically anyway), the source could change
    if (dependentData.equals(sourceData)) {
      newData.putAll(dependentData);
      return new PolicyInput(dependent.getPolicyName(), newData);
    }

    return null;
  }

  /**
   * Removing unknowns uses an even simpler algorithm: always return null, meaning the policy is
   * deleted.
   *
   * @param target existing policy
   * @param removePolicy policy to remove
   * @return null
   */
  @Override
  public PolicyInput remove(PolicyInput target, PolicyInput removePolicy) {
    return null;
  }
}