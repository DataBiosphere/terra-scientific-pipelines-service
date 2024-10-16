version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow UpdateVcfDictionaryHeader {
    input {
        File input_vcf
        File input_vcf_index
        File ref_dict
    }

    call UpdateVcfDictionaryHeader {
        input:
            vcf = input_vcf,
            vcf_index = input_vcf_index,
            ref_dict = ref_dict,
            basename = basename(input_vcf, '.vcf.gz')
    }

    output {
        File output_vcf = UpdateVcfDictionaryHeader.output_vcf
        File output_vcf_index = UpdateVcfDictionaryHeader.output_vcf_index
    }
}

task UpdateVcfDictionaryHeader {
    input {
        File vcf
        File vcf_index
        File ref_dict
        String basename

        Int disk_size_gb = ceil(4*(size(vcf, "GiB") + size(vcf_index, "GiB"))) + 20
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
        Int cpu = 1
        Int memory_mb = 6000
        Int preemptible = 0
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    command <<<
        set -e -o pipefail

        ln -sf ~{vcf} input.vcf.gz
        ln -sf ~{vcf_index} input.vcf.gz.tbi

        bcftools view -h --no-version input.vcf.gz > old_header.vcf
#        java -jar /picard.jar UpdateVcfSequenceDictionary -I old_header.vcf --SD ~{ref_dict} -O new_header.vcf

        ## update the header
        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        UpdateVCFSequenceDictionary \
        -O new_header.vcf \
        --source-dictionary ~{ref_dict} \
        --replace -V old_header.vcf \
        --disable-sequence-dictionary-validation

        bcftools reheader -h new_header.vcf -o ~{basename}.vcf.gz input.vcf.gz
        tabix ~{basename}.vcf.gz
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
