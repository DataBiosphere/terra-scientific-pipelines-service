# WDLs used to generate the reference panel files for low pass WGS / BGE imputation using GLIMPSE2

## Glimpse2SplitReference

### Purpose
This wdl takes an input reference panel vcf and splits it into separate binary 
files by chunks of chromosome. These files are used as input to GLIMPSE2 for 
imputation. This wdl is used for the reference panel creation for low pass 
WGS / BGE imputation.

#### Inputs
* contig_name
* contig_reference_chunks
* reference_filename
* reference_filename_index
* genetic_map_path_prefix
* genetic_map_path_suffix

#### Outputs
* reference_chunks
  * chunked reference `.bin` files representing the format accepted by glimpse

## FixSitesVcf

### Purpose
This wdl takes the reference panel vcf (or ideally a sites-only vcf) and performs
the following operations on it:
1) Drops all sample columns, creating a sites-only vcf
2) Trims alleles at sites with more than 30 ALTs to the most frequent 30 ALTs
3) Drops allele-specific INFO fields
4) Creates a sites table for bcftools call to use

#### Inputs
* input_vcf
* output_filename

#### Outputs
* output_vcf
* output_vcf_index
* output_table
* output_table_index

## CopyReferenceFilesToSameGcsDirectory

### Purpose
This wdl copies all the reference panel files to the same GCS directory to be used 
by the Glimpse2LowPassImputation wdl.

#### Inputs
* contigs
* sites_vcfs
* sites_vcf_indices
* sites_tables
* sites_table_indices
* reference_chunks
* google_cloud_path

#### Outputs
* SplitReferenceChunksAndCopy_stdout
* SplitReferenceChunksAndCopy_stderr
* RenameReferenceFilesAndCopy_stdout
* RenameReferenceFilesAndCopy_stderr
