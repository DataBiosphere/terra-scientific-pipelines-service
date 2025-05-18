version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow SplitMultiallelics {
    input {
        File input_vcf
        File input_vcf_index
    }

    call SeparateMultiallelics {
        input:
            original_vcf = input_vcf,
            original_vcf_index = input_vcf_index,
    }

    call CreateVcfIndex {
        input:
            vcf_input = SeparateMultiallelics.output_vcf
    }

    output {
        File multi_allelics_split_vcf = CreateVcfIndex.output_vcf
        File multi_allelics_split_vcf_index = CreateVcfIndex.output_vcf_index
    }
}

task SeparateMultiallelics {
    input {
        File original_vcf
        File original_vcf_index

        Int disk_size_gb =  ceil(2.5 * (size(original_vcf, "GiB") + size(original_vcf_index, "GiB"))) + 10
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int cpu = 1
        Int memory_mb = 6000
    }

    String output_basename = basename(original_vcf, ".vcf.gz")
    command {
        set -e -o pipefail

        bcftools norm -m - ~{original_vcf} -Oz -o ~{output_basename}.multiallelic_split.vcf.gz
    }
    output {
        File output_vcf = "~{output_basename}.multiallelic_split.vcf.gz"
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

task CreateVcfIndex {
    input {
        File vcf_input

        Int disk_size_gb = ceil(1.2 * size(vcf_input, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 6000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    }

    String vcf_basename = basename(vcf_input)

    command {
        set -e -o pipefail

        ln -sf ~{vcf_input} ~{vcf_basename}

        bcftools index -t ~{vcf_basename}
    }

    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }

    output {
        File output_vcf = "~{vcf_basename}"
        File output_vcf_index = "~{vcf_basename}.tbi"
    }
}
