version 1.0

workflow ImputationBeagleEmpty {

    String pipeline_version = "0.0.1"

    input {
        Int chunkLength = 25000000
        Int chunkOverlaps = 5000000 # this is the padding that will be added to the beginning and end of each chunk to reduce edge effects

        File multi_sample_vcf

        File ref_dict # for reheadering / adding contig lengths in the header of the ouptut VCF, and calculating contig lengths
        Array[String] contigs
        String reference_panel_path # path to the bucket where the reference panel files are stored for all contigs
        String genetic_maps_path # path to the bucket where genetic maps are stored for all contigs
        String output_callset_name # the output callset name
        Boolean split_output_to_single_sample = false

        # file extensions used to find reference panel files
        String interval_list_suffix = ".interval_list"
        String bref3_suffix = ".bref3"
    }

    call WriteEmptyFile {}

    output {
        File imputed_multi_sample_vcf = WriteEmptyFile.empty_file
        File imputed_multi_sample_vcf_index = WriteEmptyFile.empty_file
        File chunks_info = WriteEmptyFile.empty_file
        File failed_chunks = WriteEmptyFile.empty_file
        Int n_failed_chunks = 0
    }
}

task WriteEmptyFile {
    String ubuntu_docker = "ubuntu:20.04"

    command {
        touch empty_file
    }

    runtime {
        docker: ubuntu_docker
        disk: "10 GB"
        memory: "1000 MiB"
        cpu: 1
    }
    output {
        File empty_file = "empty_file"
    }
}