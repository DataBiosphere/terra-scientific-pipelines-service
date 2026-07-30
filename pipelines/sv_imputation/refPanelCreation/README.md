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

Each workflow in this directory accepts an optional `copy_to_cloud_dest` input: a
`gs://` path (bucket and prefix) that, when provided, causes the workflow's outputs to
be copied to that location via the shared `CopyToCloud` task (`CopyToCloud.wdl`). In that
case the workflow's outputs resolve to the final `gs://` destination of the files rather
than the local execution paths.

## CreateBubbleIdVcf
### Purpose
This wdl extracts a unique list of all INFO/ID values from a biallelic, sites-only bcf,
then uses that list to filter a separate input panel bcf down to only those records,
producing a vcf.

#### Inputs
* biallelic_sites_only_bcf
* biallelic_sites_only_bcf_index
* input_panel_id_split_vcf - the bcf to filter by the extracted ids
* input_panel_id_split_vcf_index
* contig - what chromosome to process of the bcf
* output_basename - base name of the final vcf
* copy_to_cloud_dest - optional gs:// path to copy the outputs to

#### Outputs
* output_panel_id_split_vcf
* output_panel_id_split_vcf_index

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
* copy_to_cloud_dest - optional gs:// path to copy the outputs to

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
* copy_to_cloud_dest - optional gs:// path to copy the outputs to

#### Outputs
* sites_only_bcf
* sites_only_bcf_index

## RelocateJsonFiles
### Purpose
This wdl takes an arbitrary JSON manifest containing `gs://` file paths (nested at any
depth, inside objects and/or arrays), copies or moves each referenced file to a single
destination `gs://` directory, and rewrites the manifest to point at the new locations.
Files are relocated by basename only (`destination_gcs_path/<basename>`); if two distinct
source paths would collide on the same destination basename, the workflow fails rather
than silently overwriting one of them.

#### Inputs
* input_json - the JSON manifest containing the `gs://` file paths to relocate
* destination_gcs_path - the `gs://` directory to relocate files to
* move_files - if true, moves (deletes the source) instead of copying; defaults to false
* dry_run - if true, computes the updated manifest and FOFN but does not actually copy or
  move any files; defaults to false

#### Outputs
* updated_json - a copy of the input manifest with all relocated paths updated to their
  new destination
* original_paths_fofn - a file of file names (FOFN) listing every original `gs://` data
  path found in the input manifest, one per line, before any relocation

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
* copy_to_cloud_dest - optional gs:// path to copy the outputs to

#### Outputs
* multi_allelics_split_bcf
* multi_allelics_split_bcf_index
