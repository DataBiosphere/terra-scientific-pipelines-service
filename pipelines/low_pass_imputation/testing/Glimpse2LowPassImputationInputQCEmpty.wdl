version 1.0

workflow InputQC {
    String pipeline_version = "0.0.1"

    input {
        # user provided inputs
        String output_basename
        File cram_manifest

        # service provided inputs
        Array[String] contigs
        String reference_panel_prefix
        File fasta
        File fasta_index
        File ref_dict

        String? pipeline_header_line # optional additional header lines to add to the output VCF
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
