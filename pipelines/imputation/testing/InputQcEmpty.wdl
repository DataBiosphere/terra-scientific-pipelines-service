version 1.0

workflow InputQC {
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

        # file extensions used to find reference panel files
        String interval_list_suffix = ".interval_list"
        String bref3_suffix = ".bref3"
    }

    call ReturnBoolAndString

    output {
        Boolean passes_qc = ReturnBoolAndString.passes_qc
        String qc_messages = ReturnBoolAndString.qc_messages
    }
}


task ReturnBoolAndString {
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
        Boolean passes_qc = true
        String qc_messages = ""
    }
}
