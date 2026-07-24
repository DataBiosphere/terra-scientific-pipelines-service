version 1.0

import "CopyToCloud.wdl" as CopyToCloudTask

# This script is under review. It is not actively tested or maintained at this time.
workflow MakeSitesOnly {
    input {
        File input_bcf
        File input_bcf_index
        String contig
        String output_basename
        String? post_contig_string
        String? copy_to_cloud_dest
    }

    call DropGenotypes {
        input:
            input_bcf = input_bcf,
            input_bcf_index = input_bcf_index,
            output_basename = "~{output_basename}.~{contig}~{post_contig_string}.sites_only"
    }

    if (defined(copy_to_cloud_dest)) {
        call CopyToCloudTask.CopyToCloud as CopyToCloud {
            input:
                source_file = DropGenotypes.output_bcf,
                source_file_index = DropGenotypes.output_bcf_index,
                copy_to_cloud_dest = select_first([copy_to_cloud_dest])
        }
    }

    output {
        String sites_only_bcf = select_first([CopyToCloud.copied_file, DropGenotypes.output_bcf])
        String sites_only_bcf_index = select_first([CopyToCloud.copied_file_index, DropGenotypes.output_bcf_index])
    }
}

task DropGenotypes {
    input {
        File input_bcf
        File input_bcf_index
        String output_basename

        Int disk_size_gb = ceil(2 * (size(input_bcf, "GiB") + size(input_bcf_index, "GiB"))) + 20
        Int cpu = 1
        Int memory_mb = 12000
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
    }

    command <<<
        set -euo pipefail

        # drop genotype (sample) columns, keeping sites only
        bcftools view -G ~{input_bcf} -Ob -o ~{output_basename}.bcf

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
        File output_bcf = "~{output_basename}.bcf"
        File output_bcf_index = "~{output_basename}.bcf.csi"
    }
}
