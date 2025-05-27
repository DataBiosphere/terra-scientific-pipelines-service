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

## BeagleImputationValidation
### Docker
File - [BeagleImputationValidationDockerFile](BeagleImputationValidationDockerFile)
Be sure to update the jars in `beagle_jars` if you
want to create a new image with updated jars

When building the docker image and you are on a M1 mac machine
, be sure to add `--platform=linux/amd64` to your docker build command

example command
```
docker build --platform=linux/amd64 -t {docker_repo}/{docker_image_name}:{docker_tag} -f {path/to/BeagleImputationValidationDockerFile} {path/to/scientificValidation/directory}
```

### Purpose
This wdl is meant to help scientifically validate
an imputation output.  This outputs a tsv file that
contains data that can be used to see how well imputation
performed like r^2 and binned by AC.

#### Inputs
* ref_panel_vcf
* ref_panel_vcf_index
* truth_vcf - vcf containing wgs data for sample(s)
that was imputed 
* truth_vcf_index
* test_vcf - output from imputation workflow
* test_vcf_index

#### Outputs
* imputed_r2_output - tsv with the following description
```
Six fields are printed for each frequency bin:

 1) MIN_AC:  the minimum reference allele count for this frequency bin
 2) MAX_AC:  the maximum reference allele count for this frequency bin
 3) ALLELES: the number of ALT alleles in this frequency bin
 4) ALT_GT:  the number of imputed genotypes for which the true genotype
             carries an ALT allele in this frequency bin
 5) DISCORD: the ALT allele dose discordance rate between imputed and true
             genotypes carrying an ALT allele in this frequency bin
 6) R2:      the squared correlation between imputed and true ALT allele dose
             for genotypes carrying an ALT alleles in this frequency bin
 ```

## LiftoverVcfs
### Purpose
This wdl takes an input vcf and lifts it over to a new reference using the gatk LiftoverVcfs tool. It is intended to be used to liftover hg19 to hg38.

#### Inputs
* vcf_path
* vcf_index_path
* liftover_chain
* hg38_reference_fasta
* hg38_reference_fasta_index
* hg38_reference_dict

Note that max_retries and preemtible_tries are set to 0 each but can be set by the user. 

#### Outputs
* hg38_vcf
* hg38_vcf_index

## SplitMultiallelics
### Purpose
This wdl takes an input vcf splits all multialleic sites into multiple biallelic records.  This is used to generate
a reference panel vcf that doesnt contain multallelic sites as Beagle doesnt know how to match a snp in the input to
a multiallelic site in the reference panel that contains that snp.

#### Inputs
* vcf_path
* vcf_index_path
* basename

Note that preemtible_tries are set to 0 each but can be set by the user.

#### Outputs
* only_biallelic_records_vcf
* only_biallelic_records_vcf_index
