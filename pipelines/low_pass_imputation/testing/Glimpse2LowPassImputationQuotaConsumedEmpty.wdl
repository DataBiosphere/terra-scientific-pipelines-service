version 1.0

workflow QuotaConsumed {
    String pipeline_version = "0.0.1"

    input {
        # user provided inputs
        String output_basename
        File cram_manifest
        Float? info_filter_for_inclusion

        # service provided inputs
        Array[String] contigs
        String reference_panel_prefix
        File fasta
        File fasta_index
        File ref_dict
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
