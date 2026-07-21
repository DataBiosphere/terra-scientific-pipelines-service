# WDLs used to generate the reference panel files for SV imputation

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
* output_basename - base name of the final vcf

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
