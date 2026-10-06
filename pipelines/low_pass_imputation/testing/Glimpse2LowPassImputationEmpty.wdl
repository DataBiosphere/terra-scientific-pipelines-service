version 1.0

workflow Glimpse2LowPassImputation {
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

        # optional additional header line to add to the output VCF
        String? pipeline_header_line
    }

    call WriteEmptyFileArray

    output {
        Array[File] imputed_vcfs = WriteEmptyFileArray.empty_files
        Array[File] imputed_vcf_indexes = WriteEmptyFileArray.empty_files
        Array[File] imputed_vcf_md5sums = WriteEmptyFileArray.empty_files

        File qc_metrics = WriteEmptyFileArray.empty_files[0]
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
