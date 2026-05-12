## StdPopSim
### Docker Image
File - [StdPopSimDockerfile](StdPopSimDockerfile)

When building the docker image and you are on a M1 mac machine, be sure to add
`--platform=linux/amd64` to your docker build command

Example command
```
docker build --platform=linux/amd64 -t {docker_repo}/{docker_image_name}:{docker_tag} -f {path/to/StdPopSimDockerFile} .
```

### Purpose

This wdl generates a simulated vcf with user provided numbers of AFR and EUR
samples using the HapMapII_GRCh38 demographic model

#### Inputs
* numAfricanSamples
* numEuropeanSamples 
* basename 
* contigs
* disk_size_gb_make_tree 
* disk_size_gb_make_vcf

#### Outputs
* vcf
* vcf_indices


## DuplicateVcfAndMerge

### Purpose
This wdl is meant to help with generating a multisample
vcf input for the imputation pipeline.  It takes in
a single vcf and a duplicate number and generates
a multisample vcf that has that number of copies of the
input vcf.  This idea for this is to use this wdl
to create the initial "small" multisample vcf that then
gets combined manually after.  This is an old version of
bcftools and doesnt perform well. YMMV

#### Inputs
* vcf_file - input vcf that will be duplicated and merged
* vcf_file_index - index for vcf_file
* vcf_duplicate_count - number of how many times to duplicate input vcf

#### Outputs
* merged_vcf - duplicated and merged vcf
* merged_vcf_index - duplicated and merged vcf index
