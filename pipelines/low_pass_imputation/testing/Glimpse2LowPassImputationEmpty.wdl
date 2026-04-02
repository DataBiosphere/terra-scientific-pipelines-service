version 1.0

workflow Glimpse2LowPassImputation {
    input {
        String pipeline_version = "0.0.1"

        Array[String] contigs

        String reference_panel_prefix

        File? input_vcf
        File? input_vcf_index
        Array[File]? crams
        Array[File]? cram_indices
        Array[String] sample_ids
        File fasta
        File fasta_index
        String output_basename
        File? cram_manifest

        File ref_dict

        Boolean impute_reference_only_variants = false
        Boolean call_indels = false

        Int calling_batch_size = 100

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.0.0"
        String glimpse_docker = "us.gcr.io/broad-dsde-methods/glimpse:kachulis_ck_bam_reader_retry_cf5822c"
    }

    call WriteEmptyFile

    output {
        File imputed_vcf = WriteEmptyFile.empty_file
        File imputed_vcf_index = WriteEmptyFile.empty_file
        File imputed_vcf_md5sum = WriteEmptyFile.empty_file

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
