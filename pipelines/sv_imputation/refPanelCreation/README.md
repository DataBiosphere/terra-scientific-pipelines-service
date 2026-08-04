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

## CreateRefPanelSampleSubset
### Purpose
This wdl subsets a parent SV reference panel bcf down to a given sample list for a single
chromosome, then derives the companion files needed downstream: it splits multiallelic
sites into biallelic records (recalculating AC/AN and dropping AC=0 sites), drops
genotype columns to produce a sites-only bcf, and filters a companion id-split vcf down
to only the INFO/ID values present in that sites-only bcf.

#### Inputs
* input_panel_bubble_bcf - source (bubble-annotated) reference panel bcf to subset
* input_panel_bubble_bcf_index - index for the above input
* input_panel_id_split_vcf - companion id-split vcf to filter down to the ids retained in
  the generated sites-only bcf
* input_panel_id_split_vcf_index - index for the above input
* sample_list - list of samples to subset input_panel_bubble_bcf to
* contig - what chromosome is being processed
* output_basename - base name used for all generated output files

#### Outputs
* bubble_split_bcf - sample-subset bcf with multiallelic sites split to biallelic
* bubble_split_bcf_index - index of the above
* bubble_split_sites_bcf - sites-only (genotype-dropped) version of bubble_split_bcf
* bubble_split_sites_bcf_index - index of the above
* panel_id_split_vcf - input_panel_id_split_vcf filtered down to the ids present in
  bubble_split_sites_bcf
* panel_id_split_vcf_index - index of the above

## RelocateAllSVReferencePanelFiles
### Purpose
This wdl relocates all of the final SV reference panel outputs to their destinations:
every `gs://` file referenced by the `chunked_panel_json` and `pop_panel_resources_json`
manifests (via the shared `RelocateFiles` task from `RelocateJsonFiles.wdl`), both
manifests themselves, and the `preprocess_panel_bubble_split_sites_only_vcf` file + index
produced separately by `CreatePanelAuxiliaryFiles`. `preprocess_panel_bubble_split_sites_only_vcf`
and its index are passed as `gs://` path Strings (not `File`) so that a move actually
deletes the source object, rather than only deleting a Cromwell-localized local copy.

#### Inputs
* chunked_panel_json - manifest containing the `gs://` paths of the chunked panel bin
  files to relocate
* pop_panel_resources_json - the `pop_panel_resources_json` manifest containing the
  `gs://` file paths to relocate
* preprocess_panel_bubble_split_sites_only_vcf - the sites-only vcf produced by
  CreatePanelAuxiliaryFiles, as a `gs://` path
* preprocess_panel_bubble_split_sites_only_vcf_idx - index for the above input, as a
  `gs://` path
* panel_bin_files_destination_gcs_path - the `gs://` directory to relocate the chunked
  panel bin files to
* auxiliary_files_destination_gcs_path - the `gs://` directory to relocate the
  pop_panel_resources files and the sites-only vcf + index to
* json_destination_gcs_path - the `gs://` directory both rewritten manifests are copied to
* move_files - if true, moves (deletes the source) instead of copying; defaults to false
* dry_run - if true, computes the updated manifests and FOFNs but does not actually copy
  or move any files; defaults to false

#### Outputs
* chunked_panel_updated_json / pop_panel_resources_updated_json - local copies of the
  rewritten manifests with all relocated paths updated to their new destination
* chunked_panel_original_paths_fofn / pop_panel_resources_original_paths_fofn - a file of
  file names (FOFN) listing every original `gs://` data path found in each manifest, one
  per line, before any relocation
* chunked_panel_updated_json_gcs_path / pop_panel_resources_updated_json_gcs_path - the
  `gs://` path each rewritten manifest was copied to
* relocated_preprocess_panel_bubble_split_sites_only_vcf - the `gs://` path the sites-only
  vcf was relocated to
* relocated_preprocess_panel_bubble_split_sites_only_vcf_idx - the `gs://` path the index
  was relocated to

# Support WDLs whose tasks are used in above WDLs but can also be run directly as needed

## MakeSitesOnly
### Purpose
This wdl takes an input bcf and drops all sample information, producing a sites-only bcf.

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
than silently overwriting one of them. The rewritten manifest itself (under its original
basename) is then copied to `json_destination_gcs_path`.

#### Inputs
* input_json - the JSON manifest containing the `gs://` file paths to relocate
* destination_gcs_path - the `gs://` directory to relocate files to
* json_destination_gcs_path - the `gs://` directory the rewritten manifest is copied to
* move_files - if true, moves (deletes the source) instead of copying; defaults to false
* dry_run - if true, computes the updated manifest and FOFN but does not actually copy or
  move any files, including the manifest itself; defaults to false

#### Outputs
* updated_json - a local copy of the rewritten manifest with all relocated paths updated
  to their new destination
* original_paths_fofn - a file of file names (FOFN) listing every original `gs://` data
  path found in the input manifest, one per line, before any relocation
* updated_json_gcs_path - the `gs://` path the rewritten manifest was copied to

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

