
## CreateImputationRefPanelBeagle

### Purpose
This wdl takes in ref panel vcf for a chromosome and creates bref files,
interval lists, bed files, and unique variant id lists from that vcf.

#### Inputs
* ref_vcf - ref panel vcf for a chromosome
* ref_vcf_index - corresponding ref panel vcf index
* ref_dict - reference dictionary for the reference
* create_brefs - boolean to generate bref files, defaults to true
* create_interval_lists - boolean to generate interval list files, defaults to true
* create_bed_files - boolean to generate bed files, defaults to true
* create_unique_variant_id_lists - boolean to generate unique variant id lists, defaults to true

#### Outputs
* interval_list
* bed_file
* bref
* unique_variant_ids


## ReshapeReferencePanel
### Purpose
This wdl applies the reidentifying mitigation algorithm
to a reference panel vcf file.  The logic is described
in this repo - https://github.com/TheoCavinato/RESHAPE.
This wdl parallelizes over intervals and samples in
order to speed up the wallclock time of the workflow.

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
* output_vcf - output vcf after mitigation algorithm has been run
* output_vcf_index - index of output vcf

## SplitMultiallelics
### Purpose
This wdl takes an input vcf and splits all multialleic sites into multiple biallelic records.  This is used to generate
a reference panel vcf that doesn't contain multiallelic sites as Beagle doesnt know how to match a snp in the input to
a multiallelic site in the reference panel that contains that snp.

#### Inputs
* vcf_path
* vcf_index_path
* basename

Note that preemptible_tries are set to 0 each but can be set by the user.

#### Outputs
* only_biallelic_records_vcf
* only_biallelic_records_vcf_index

