version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow SplitMultiallelics {
    input {
        File input_vcf
        File input_vcf_index
        String basename
    }

    call SeparateMultiallelics {
        input:
            original_vcf = input_vcf,
            original_vcf_index = input_vcf_index,
            output_basename = basename + ".only_biallelic_records"
    }

    output {
        File only_biallelic_records_vcf = SeparateMultiallelics.output_vcf
        File only_biallelic_records_vcf_index = SeparateMultiallelics.output_vcf_index
    }
}

task SeparateMultiallelics {
    input {
        File original_vcf
        File original_vcf_index
        String output_basename

        Int disk_size_gb =  ceil(2*(size(original_vcf, "GiB") + size(original_vcf_index, "GiB"))) + 10
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int cpu = 1
        Int memory_mb = 6000
    }
    command {
        set -e -o pipefail

        bcftools norm -m - ~{original_vcf} -Oz -o ~{output_basename}.vcf.gz
        bcftools index -t ~{output_basename}.vcf.gz
    }
    output {
        File output_vcf = "~{output_basename}.vcf.gz"
        File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
    }
    runtime {
        docker: bcftools_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 0
        noAddress: true
    }
}
