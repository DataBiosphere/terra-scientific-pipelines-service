version 1.0

# name it ImputationBeagle even though it's empty for testing
workflow ImputationBeagle {

    String pipeline_version = "0.0.1"

    input {
        Int chunkLength = 25000000
        Int chunkOverlaps = 5000000

        File multi_sample_vcf

        File ref_dict
        Array[String] contigs
        String reference_panel_path_prefix
        String genetic_maps_path
        String output_basename

        String? pipeline_header_line
        Float? min_dr2_for_inclusion

        # file extensions used to find reference panel files
        String interval_list_suffix = ".interval_list"
        String bref3_suffix = ".bref3"
    }

    call WriteEmptyFile

    output {
        File imputed_multi_sample_vcf = WriteEmptyFile.empty_file
        File imputed_multi_sample_vcf_index = WriteEmptyFile.empty_file
        File chunks_info = WriteEmptyFile.empty_file
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
        maxRetries: 2
    }
    output {
        File empty_file = "empty_file"
    }
}
