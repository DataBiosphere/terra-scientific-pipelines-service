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
        File subset_vcf = UpdateHeader.output_vcf
        File subset_vcf_index = UpdateHeader.output_vcf_index
    }
}

task UpdateHeader {
    input {
        File vcf
        File vcf_index
        File ref_dict
        String basename

        Int disk_size_gb = ceil(4*(size(vcf, "GiB") + size(vcf_index, "GiB"))) + 20
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
        Int cpu = 1
        Int memory_mb = 6000
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    command <<<

        ## update the header of the merged vcf
        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        UpdateVCFSequenceDictionary \
        --source-dictionary ~{ref_dict} \
        --output ~{basename}.vcf.gz \
        --replace -V ~{vcf} \
        --disable-sequence-dictionary-validation
    >>>
    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 3
    }
    output {
        File output_vcf = "~{basename}.vcf.gz"
        File output_vcf_index = "~{basename}.vcf.gz.tbi"
    }
}
