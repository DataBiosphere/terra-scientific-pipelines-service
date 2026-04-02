version 1.0

workflow InputQC {
    String pipeline_version = "0.0.1"

    input {
        Array[String] contigs

        String reference_panel_prefix

        Array[File]? crams
        Array[File]? cram_indices
        Array[String]? sample_ids
        File? cram_manifest
        File fasta
        File fasta_index
        String output_basename

        File ref_dict
    }

    call MockValidateInputs

    output {
        Boolean passes_qc = MockValidateInputs.passes_qc
        String qc_messages = MockValidateInputs.qc_messages
    }
}

task MockValidateInputs  {
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