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


## SubsetVcfByBedFile
### Purpose
This wdl is meant to be used to subset a vcf down
to sites provided through a bed file.  This wdl does
not interact with headers or annotations mostly because
the only really "required" header is the dictionary
and that gets transferred across and the imputation
tool only look at GT and no info/format fields so
we can just leave them be.

This wdl can optionally additionally extract a subset of samples if the
optional input `samples_to_subselect` is provided.

#### Inputs
* input_vcf - input file to be subset
* input_vcf_index 
* bed_file - bed file containing intervals to subset by
* (optional) samples_to_subselect - text file containing 
newline-separated list of samples to subset the vcf by

#### Outputs
* subset_vcf - subsetted vcf
* subset_vcf_index 

## UpdateVcfDictionaryHeader
### Purpose
This wdl can be used to update the dictionary header lines
(i.e. `##contigs=` lines) in a vcf file, based on a supplied
dictionary file. This is useful when the contigs in the vcf
file do not match the contigs in the reference panel used
for imputation. The optional input `disable_sequence_dictionary_validation`,
by default set to true, checks the contigs in the input vcf file and
throws an error if they are not valid. If set to false, the contigs
in the input vcf file are not checked.

#### Inputs
* input_vcf - input vcf file
* input_vcf_index
* dictionary_file - dictionary file to use for updating the vcf header;
we usually use `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict`.

#### Outputs
* output_vcf - reheadered vcf file
* output_vcf_index
