## ReshapeReferencePanel
### Purpose
This wdl applies the reidentifying mitigation algorithm
to a reference panel vcf file.  The logic is described
in this repo - https://github.com/TheoCavinato/RESHAPE.
This wdl is basically a wrapper around that tool/image

#### Inputs
* ref_panel_vcf - input reference panel vcf to mitigate
* genetic_map - recombination map
* chromosome - what chromosome process of vcf

#### Outputs
* recombined_reference_panel - output vcf after mitigation 
algorithm has been run
