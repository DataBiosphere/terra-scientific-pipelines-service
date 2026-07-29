version 1.0

import "CopyToCloud.wdl" as CopyToCloudTask

# This script is under review. It is not actively tested or maintained at this time.
workflow CreateSimpleIdBcf {
    input {
        String contig
        File biallelic_sites_only_bcf
        File biallelic_sites_only_bcf_index
        File input_panel_id_split_vcf
        File input_panel_id_split_vcf_index
        String output_basename
        String? copy_to_cloud_dest
    }

    call ExtractIdsAndFilter {
        input:
            biallelic_sites_only_bcf = biallelic_sites_only_bcf,
            biallelic_sites_only_bcf_index = biallelic_sites_only_bcf_index,
            input_panel_id_split_vcf = input_panel_id_split_vcf,
            input_panel_id_split_vcf_index = input_panel_id_split_vcf_index,
            output_basename = "~{output_basename}.~{contig}.split.sites.simple"
    }

    if (defined(copy_to_cloud_dest)) {
        call CopyToCloudTask.CopyToCloud as CopyToCloud {
            input:
                source_file = ExtractIdsAndFilter.output_panel_simple_id_split_bcf,
                source_file_index = ExtractIdsAndFilter.output_panel_simple_id_split_bcf_index,
                copy_to_cloud_dest = select_first([copy_to_cloud_dest])
        }
    }

    output {
        String output_panel_simple_id_split_bcf = select_first([CopyToCloud.copied_file, ExtractIdsAndFilter.output_panel_simple_id_split_bcf])
        String output_panel_simple_id_split_bcf_index = select_first([CopyToCloud.copied_file_index, ExtractIdsAndFilter.output_panel_simple_id_split_bcf_index])
    }
}


task ExtractIdsAndFilter {
    input {
        File biallelic_sites_only_bcf
        File biallelic_sites_only_bcf_index
        File input_panel_id_split_vcf
        File input_panel_id_split_vcf_index
        String output_basename

        Int disk_size_gb = ceil(2 * (size(input_panel_id_split_vcf, "GiB") + size(input_panel_id_split_vcf_index, "GiB"))) + 20
        Int cpu = 1
        Int memory_mb = 12000
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
    }


    command <<<
        set -euo pipefail

        echo "doing query"
        bcftools query -f '%INFO/ID\n' ~{biallelic_sites_only_bcf} | grep -v : > simple.ids.list

        echo "doing view"
        bcftools view -i "INFO/ID=@simple.ids.list" -Ob -o ~{output_basename}.bcf ~{input_panel_id_split_vcf}

        bcftools index ~{output_basename}.bcf
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
        File output_panel_simple_id_split_bcf = "~{output_basename}.bcf"
        File output_panel_simple_id_split_bcf_index = "~{output_basename}.bcf.csi"
    }
}
