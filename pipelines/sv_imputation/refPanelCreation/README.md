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

## GeneratePreprocessPanelBubbleSplitSitesOnlyBcf
### Purpose
This wdl is used to create the contig specific
`preprocess_panel_bubble_split_sites_only_vcf` input to the
SV imputation wdl from the contig specific `reduced_panel_bubble_vcf`
file produced from the PostprocessBubblePanel wdl.

#### Inputs
* reduced_panel_bubble_vcf - reduced_panel_bubble_vcf output from PostprocessBubblePanel wdl
* reduced_panel_bubble_vcf_index - index file
* contig - what chromosome is being processed
* output_base_name - base name of final reduced sites only bcf file

#### Outputs
* output_bcf - output bcf after 
* output_bcf_index - index of output bcf

## CreatePanelAuxiliaryFiles
### Purpose
This wdl is used to create the `pop_glimpse2_panel_resources_json`
and `preprocess_panel_bubble_split_sites_only_vcf` + index file used
as input to sv imputation.  The inputs to these files are outputs
produced from PostprocessBubblePanel and GeneratePreprocessPanelBubbleSplitSitesOnlyBcf

#### Inputs
* chromosomes - chromosomes to process, this will be the first key in each json object
* panel_bubble_split_sites_only_vcf_array - supplied from the `panel_bubble_split_sites_only_vcf` output of PostprocessBubblePanel
* panel_bubble_split_sites_only_vcf_idx_array - index for the above input
* panel_id_split_vcf_gz_array - supplied from the `panel_id_split_vcf_gz` output of PostprocessBubblePanel
* panel_id_split_vcf_gz_tbi_array - index for the above input
* reduced_panel_bubble_split_simple_sites_bcf - supplied from the `output_bcf` output of GeneratePreprocessPanelBubbleSplitSitesOnlyBcf
* reduced_panel_bubble_split_simple_sites_bcf_index - index for the above input
* output_prefix - base name of the `gathered_reduced_panel_bubble_split_simple_sites_bcf` output
* ref_dict - reference dictionary to get chromosome length from so get proper pop_regions found in output json
* interval_size - interval size for pop_regions found in output json

#### Outputs
* pop_panel_resources_json - json used as `pop_panel_resources_json` to sv imputation workflow
* gathered_reduced_panel_bubble_split_simple_sites_bcf - file used as `preprocess_panel_bubble_split_sites_only_vcf` input to sv imputation workflow
* gathered_reduced_panel_bubble_split_simple_sites_bcf_index - file used as `preprocess_panel_bubble_split_sites_only_vcf_idx` input to sv imputation workflow


# WDLs used to generate the testing reference panel files for SV imputation

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

## SplitMultiallelicsBcf
### Purpose
This wdl takes an input bcf and splits all multiallelic sites into multiple biallelic
records. It parallelizes over base chunks of the chromosome to speed up wallclock time.
This is the bcf equivalent of the array imputation SplitMultiallelics wdl.

#### Inputs
* input_bcf
* input_bcf_index
* ref_dict - reference dictionary for the reference
* contig - what chromosome to process of the bcf
* output_basename - base name of the final bcf
* num_base_chunk_size - number of bases to split the chromosome by (default of 10000000)

#### Outputs
* multi_allelics_split_bcf
* multi_allelics_split_bcf_index
