# WDLs used to generate the reference panel files for SV imputation using GLIMPSE2

## ReshapeReferencePanel
### Purpose
This wdl applies the reidentifying mitigation algorithm
to a reference panel bcf file.  The logic is described
in this repo - https://github.com/TheoCavinato/RESHAPE.

#### Inputs
* ref_panel_vcf - input reference panel vcf to mitigate
* genetic_map - recombination map
* contig - what chromosome to process of vcf
* output_base_name - base name of final reshaped vcf
* ref_dict - reference dictionary file
* reshape_threads - number of threads to use for reshape tool
* num_base_chunk_size - number of bases to split the chromosome by (default of 25000000)
* sample_chunk_size - number of samples to split by (default of 50000)

#### Outputs
* output_bcf - output bcf after mitigation algorithm has been run
* output_bcf_index - index of output bcf
