version 1.0

import "CopyToCloud.wdl" as CopyToCloudTask

# This script is under review. It is not actively tested or maintained at this time.
workflow CreateBubbleIdVcf {
    input {
        String contig
        File biallelic_sites_only_bcf
        File biallelic_sites_only_bcf_index
        File input_panel_id_split_vcf
        File input_panel_id_split_vcf_index
        String output_basename
        String? copy_to_cloud_dest
    }

    call ExtractIds {
        input:
            biallelic_sites_only_bcf = biallelic_sites_only_bcf,
            biallelic_sites_only_bcf_index = biallelic_sites_only_bcf_index
    }

    call FilterByIds {
        input:
            input_panel_id_split_vcf = input_panel_id_split_vcf,
            input_panel_id_split_vcf_index = input_panel_id_split_vcf_index,
            ids_to_keep = ExtractIds.ids_to_keep,
            output_basename = "~{output_basename}.~{contig}.id.split"
    }

    if (defined(copy_to_cloud_dest)) {
        call CopyToCloudTask.CopyToCloud as CopyToCloud {
            input:
                source_file = FilterByIds.output_panel_id_split_vcf,
                source_file_index = FilterByIds.output_panel_id_split_vcf_index,
                copy_to_cloud_dest = select_first([copy_to_cloud_dest])
        }
    }

    output {
        String output_panel_id_split_vcf = select_first([CopyToCloud.copied_file, FilterByIds.output_panel_id_split_vcf])
        String output_panel_id_split_vcf_index = select_first([CopyToCloud.copied_file_index, FilterByIds.output_panel_id_split_vcf_index])
    }
}

task ExtractIds {
    input {
        File biallelic_sites_only_bcf
        File biallelic_sites_only_bcf_index

        Int disk_size_gb = ceil(2 * (size(biallelic_sites_only_bcf, "GiB") + size(biallelic_sites_only_bcf_index, "GiB"))) + 20
        Int cpu = 1
        Int memory_mb = 12000
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
    }

    command <<<
        set -euo pipefail

        # extract a unique list of all INFO/ID values; a single record's ID field may
        # contain multiple IDs separated by colons or commas, so split those onto
        # their own lines before deduplicating. drop empty and missing (".") entries.
        bcftools query -f '%INFO/ID\n' ~{biallelic_sites_only_bcf} \
            | tr ',:' '\n\n' \
            | sed -e '/^$/d' -e '/^\.$/d' \
            | sort -u > ids_to_keep.list
    >>>

    runtime {
        docker: bcftools_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 0
        noAddress: true
    }

    output {
        File ids_to_keep = "ids_to_keep.list"
    }
}

task FilterByIds {
    input {
        File input_panel_id_split_vcf
        File input_panel_id_split_vcf_index
        File ids_to_keep
        String output_basename

        Int disk_size_gb = ceil(2 * (size(input_panel_id_split_vcf, "GiB") + size(input_panel_id_split_vcf_index, "GiB"))) + 20
        Int cpu = 1
        Int memory_mb = 12000
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
    }


    command <<<
        set -euo pipefail

        # keep only records whose INFO/ID is in the ids_to_keep list
        bcftools view \
            -i 'INFO/ID=@~{ids_to_keep}' \
            ~{input_panel_id_split_vcf} \
            -O z -o ~{output_basename}.vcf.gz

        bcftools index -t ~{output_basename}.vcf.gz
    >>>

    runtime {
        docker: bcftools_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 0
        noAddress: true
    }

    output {
        File output_panel_id_split_vcf = "~{output_basename}.vcf.gz"
        File output_panel_id_split_vcf_index = "~{output_basename}.vcf.gz.tbi"
    }
}
