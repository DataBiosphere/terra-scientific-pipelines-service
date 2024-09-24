## ReshapeReferencePanel
### Purpose
This wdl applies the reidentifying mitigation algorithm
to a reference panel vcf file.  The logic is described
in this repo - https://github.com/TheoCavinato/RESHAPE.
This wdl is basically a wrapper around that tool/image

#### Inputs
* ref_panel_vcf - input reference panel vcf to mitigate
* genetic_map - recombination map
* chromosome - what chromosome to process of vcf

#### Outputs
* recombined_reference_panel - output vcf after mitigation 
algorithm has been run


## ReshapeReferencePanel
### Purpose
This wdl is meant to be used to subset a vcf down
to sites provided through a bed file

#### Inputs
* input_vcf - input file to be subset
* input_vcf_index 
* bed_file - bed file containing intervals to subset by

#### Outputs
* subset_vcf - subsetted vcf
* subset_vcf_index 
