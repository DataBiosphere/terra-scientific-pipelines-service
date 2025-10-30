version 1.0

# Extracts unique variant identifiers in the format %CHROM:%POS:%REF:%ALT from a VCF file

workflow ExtractUniqueVariantIds {
    input {
        File input_vcf
    }

    call ExtractUniqueVariantIds {
        input:
            vcf = input_vcf
    }

    output {
        File unique_variant_ids = ExtractUniqueVariantIds.unique_variant_ids
    }
}

task ExtractUniqueVariantIds {
    input {
        File vcf

        Int disk_size_gb = ceil(2 * size(vcf, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 4000
        String bcftools_docker = "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
        Int preemptible = 3
    }

    String output_base = basename(vcf, ".vcf.gz")

    command <<<
        set -e -o pipefail

        bcftools query -f "%CHROM:%POS:%REF:%ALT\n" ~{vcf} | uniq > ~{output_base}.unique_variants.raw
        # sort and remove blank lines
        sort ~{output_base}.unique_variants.raw | sed '/^$/d' > ~{output_base}.unique_variants
    >>>

    runtime {
        docker: bcftools_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: 1
        noAddress: true
    }

    output {
        File unique_variant_ids = "~{output_base}.unique_variants"
    }
}
