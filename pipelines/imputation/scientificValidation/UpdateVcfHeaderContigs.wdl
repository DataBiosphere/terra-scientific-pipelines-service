version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow UpdateVcfHeaderContigs {
    input {
        File input_vcf
        File input_vcf_index
        File ref_dict
    }

    call UpdateHeader {
        input:
            vcf = input_vcf,
            vcf_index = input_vcf_index,
            ref_dict = ref_dict,
            basename = basename(input_vcf, '.vcf.gz')
    }

    output {
        File output_vcf = UpdateHeader.output_vcf
        File output_vcf_index = UpdateHeader.output_vcf_index
    }
}

task UpdateHeader {
    input {
        File vcf
        File vcf_index
        File ref_dict
        String basename
        Boolean disable_sequence_dictionary_validation = true

        Int disk_size_gb = ceil(4*(size(vcf, "GiB") + size(vcf_index, "GiB"))) + 20
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
        Int cpu = 1
        Int memory_mb = 6000
        Int preemptible = 0
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    String disable_sequence_dict_validation_flag = if disable_sequence_dictionary_validation then "--disable-sequence-dictionary-validation" else ""

    command <<<
        set -e -o pipefail

        ln -sf ~{vcf} input.vcf.gz
        ln -sf ~{vcf_index} input.vcf.gz.tbi

        ## update the header of the merged vcf
        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        UpdateVCFSequenceDictionary \
        --source-dictionary ~{ref_dict} \
        --output ~{basename}.vcf.gz \
        --replace -V input.vcf.gz \
        ~{disable_sequence_dict_validation_flag}
    >>>
    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: preemptible
    }
    output {
        File output_vcf = "~{basename}.vcf.gz"
        File output_vcf_index = "~{basename}.vcf.gz.tbi"
    }
}
