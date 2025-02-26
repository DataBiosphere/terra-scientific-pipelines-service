## CompareVcfs
### Purpose
This wdl is intended to be used to compare the imputed output VCFs from different versions of the ImputationBeagle wdl. 
It runs a similar check to the [VerifyImputationBeagle.wdl in warp](https://github.com/broadinstitute/warp/blob/develop/verification/VerifyImputationBeagle.wdl). 
When used with patternForLinesToExcludeFromComparison = "##", it ignores the header lines in the VCFs but checks the sample 
names and the remaining contents of the VCFs. The workflow fails if any differences are found.

#### Inputs
* test_gvcf
* truth_gvcf

#### Outputs
none (workflow success indicates that the VCFs are identical)

## ImputationBeagleEmpty
### Purpose
This wdl is used in e2e tests for Teaspoons and its cli, and sometimes with local deployments of Teaspoons, for testing 
infrastructure without worrying about the imputation process. It takes the same inputs and returns the same outputs as the 
real [ImputationBeagle wdl in warp](https://github.com/broadinstitute/warp/tree/develop/pipelines/broad/arrays/imputation_beagle), 
but it only runs a dummy task and returns empty files.

## QuotaConsumedEmpty
### Purpose
Similar to ImputationBeagleEmpty, This wdl is used in e2e tests for Teaspoons and its cli, and sometimes with local deployments of Teaspoons, for testing
infrastructure without worrying about the quota calculation process. It takes the same inputs and returns the same outputs as the
real [ArrayImputationQuotaConsumed wdl in warp](https://github.com/broadinstitute/warp/tree/develop/pipelines/broad/arrays/imputation_beagle),
but it only runs a dummy task and returns empty files.
