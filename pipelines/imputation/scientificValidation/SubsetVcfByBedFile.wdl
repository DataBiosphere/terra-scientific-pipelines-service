version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow SubsetVcfByBedFile {
    input {
        File input_vcf
        File input_vcf_index
        File bed_file
    }

    call BcftoolsSubsetVcf {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            bed_file = bed_file
    }

    output {
        File subset_vcf = BcftoolsSubsetVcf.output_vcf
        File subset_vcf_index = BcftoolsSubsetVcf.output_vcf_index
    }
}

task BcftoolsSubsetVcf {
    input {
        File input_vcf
        File input_vcf_index
        File bed_file

        Int disk_size_gb = ceil(3 * size(input_vcf, "GiB")) + 20
        Int cpu = 1
        Int memory_mb = 6000
    }

    String basename = basename(input_vcf, '.vcf.gz')

    command {
        set -e -o pipefail

        bcftools view \
        -R ~{bed_file} \
        -O z \
        -o ~{basename}.subset.vcf.gz \
        ~{input_vcf}

        bctools index -t ~{basename}.subset.vcf.gz

    }

    output {
        File output_vcf = "~{basename}.subset.vcf.gz"
        File output_vcf_index = "~{basename}.subset.vcf.gz.tbi"
    }

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.5.0.0"
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }
}
