# WDLs used to generate the reference panel files for SV imputation

## ExtractSamplesAndFilter
### Purpose
This wdl takes an input bcf, subsets it to a given list of samples, and then
filters the result to keep only alt sites (i.e. removing hom ref sites).

#### Inputs
* input_bcf
* input_bcf_index
* sample_list - list of samples to subset the bcf to
* contig - what chromosome to process of the bcf
* output_basename - base name of the final bcf
* post_contig_string - optional string appended after the contig in the output name

#### Outputs
* output_bcf
* output_bcf_index

## MakeSitesOnly
### Purpose
This wdl takes an input bcf and drops the genotype (sample) columns, producing a
sites-only bcf.

#### Inputs
* input_bcf
* input_bcf_index
* contig - what chromosome to process of the bcf
* output_basename - base name of the final bcf
* post_contig_string - optional string appended after the contig in the output name

#### Outputs
* sites_only_bcf
* sites_only_bcf_index
