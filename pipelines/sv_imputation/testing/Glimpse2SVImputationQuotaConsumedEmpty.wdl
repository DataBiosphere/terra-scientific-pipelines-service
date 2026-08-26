version 1.0

workflow QuotaConsumed {
    input {
        # user defined inputs
        File gvcf_manifest
        String output_basename
        Float? info_filter_for_inclusion

        # service provided inputs
        File preprocess_panel_bubble_split_sites_only_vcf
        File preprocess_panel_bubble_split_sites_only_vcf_idx
        Array[String] paste_regions
        Array[String] chromosomes
        File genetic_maps_tsv
        File ref_dict
        File chunked_panel_json
        File pop_glimpse2_panel_resources_json
        String? pipeline_header_line
    }

    call CalculateMockQuota

    output {
        Int quota_consumed = CalculateMockQuota.quota_consumed
    }
}

task CalculateMockQuota  {
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
