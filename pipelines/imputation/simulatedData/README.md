## StdPopSim
### Docker Image
When building the docker image and you are on a M1 mac machine, be sure to add `--platform=linux/amd64` to your docker build command

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

## CreateImputationRefPanelBeagle

### Purpose
This wdl takes in ref panel vcfs and creates bref files,
interval lists, and bed files from those vcfs.  The vcf
input is expected to be in order by chromosome

#### Inputs
* ref_vcf - a list of ref panel vcfs in order by chromosome
* ref_vcf_index - corresponding list of ref panel vcfs
* ref_dict - reference dictionary for the reference

#### Outputs
* interval_lists
* bed_files
* brefs
