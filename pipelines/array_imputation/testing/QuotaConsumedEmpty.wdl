version 1.0

workflow QuotaConsumed {
    input {
        # user provided inputs
        File multi_sample_vcf
        String output_basename
        Float? min_dr2_for_inclusion

        # service provided inputs
        Array[String] contigs
        String genetic_maps_path
        File ref_dict
        String reference_panel_path_prefix
        String? pipeline_header_line
    }

    call ReturnHardcodedInt

    output {
        Int quota_consumed = ReturnHardcodedInt.quota_consumed
    }
}


task ReturnHardcodedInt {
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
        Int quota_consumed = 50
    }
}
