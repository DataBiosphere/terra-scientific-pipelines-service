version 1.0

workflow Glimpse2LowPassImputation {
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

    call WriteEmptyFile

    output {
        File imputed_vcf = WriteEmptyFile.empty_file
        File imputed_vcf_index = WriteEmptyFile.empty_file
        File imputed_vcf_md5sum = WriteEmptyFile.empty_file

        File imputed_hom_ref_sites_only_vcf = WriteEmptyFile.empty_file
        File imputed_hom_ref_sites_only_vcf_index = WriteEmptyFile.empty_file
        File imputed_hom_ref_sites_only_vcf_md5 = WriteEmptyFile.empty_file

        File qc_metrics = WriteEmptyFile.empty_file
        File? coverage_metrics = WriteEmptyFile.empty_file
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
