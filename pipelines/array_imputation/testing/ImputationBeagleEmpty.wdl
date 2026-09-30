version 1.0

# name it ImputationBeagle even though it's empty for testing
workflow ImputationBeagle {
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

    call WriteEmptyFileArray

    output {
        Array[File] imputed_multi_sample_vcfs = WriteEmptyFileArray.empty_files
        Array[File] imputed_multi_sample_vcf_indexes = WriteEmptyFileArray.empty_files

        File chunks_info = WriteEmptyFileArray.empty_files[0]
        File contigs_info = WriteEmptyFileArray.empty_files[0]
    }
}

task WriteEmptyFileArray {
    String ubuntu_docker = "ubuntu:20.04"

    command {
        touch empty_file_0
        touch empty_file_1
    }

    runtime {
        docker: ubuntu_docker
        disk: "10 GB"
        memory: "1000 MiB"
        cpu: 1
        maxRetries: 2
    }
    output {
        Array[File] empty_files = glob("empty_file*")
    }
}
