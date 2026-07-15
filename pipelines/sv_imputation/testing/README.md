## CompareBcfs
### Purpose
This wdl is intended to be used to compare the imputed output BCFs from different versions of the SV imputation wdl.
When used with patternForLinesToExcludeFromComparison = "##", it ignores the header lines in the BCFs but checks the sample
names and the remaining contents of the BCFs. The workflow fails if any differences are found.

#### Inputs
* test_bcf
* truth_bcf

#### Outputs
none (workflow success indicates that the BCFs are identical aside from headers)
