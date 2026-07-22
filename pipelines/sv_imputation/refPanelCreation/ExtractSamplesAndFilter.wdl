version 1.0

import "CopyToCloud.wdl" as CopyToCloudTask

# This script is under review. It is not actively tested or maintained at this time.
workflow ExtractSamplesAndFilter {
    input {
        File input_bcf
        File input_bcf_index
        File sample_list
        String contig
        String output_basename
        String? post_contig_string
        String? copy_to_cloud_dest
    }

    call ExtractAndFilter {
        input:
            input_bcf = input_bcf,
            input_bcf_index = input_bcf_index,
            sample_list = sample_list,
            contig = contig,
            output_basename = output_basename,
            post_contig_string = post_contig_string
    }

    if (defined(copy_to_cloud_dest)) {
        call CopyToCloudTask.CopyToCloud as CopyToCloud {
            input:
                source_file = ExtractAndFilter.output_bcf,
                source_file_index = ExtractAndFilter.output_bcf_index,
                copy_to_cloud_dest = select_first([copy_to_cloud_dest])
        }
    }

    output {
        String output_bcf = select_first([CopyToCloud.copied_file, ExtractAndFilter.output_bcf])
        String output_bcf_index = select_first([CopyToCloud.copied_file_index, ExtractAndFilter.output_bcf_index])
    }
}

task ExtractAndFilter {
    input {
        File input_bcf
        File input_bcf_index
        File sample_list
        String contig
        String output_basename
        String post_contig_string = ""

        Int disk_size_gb = ceil(3 * (size(input_bcf, "GiB") + size(input_bcf_index, "GiB"))) + 20
        Int cpu = 1
        Int memory_mb = 12000
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
    }

    command <<<
        set -euo pipefail

        # subset samples
        bcftools view \
            -S ~{sample_list} \
            ~{input_bcf} \
            -O b -o sample_subset.bcf

        bcftools index sample_subset.bcf

        # keep alt sites (i.e. remove hom ref sites)
        bcftools view \
            -i 'GT[*]="alt"' \
            sample_subset.bcf \
            -O b \
            -o ~{output_basename}.~{contig}~{post_contig_string}.bcf

        bcftools index ~{output_basename}.~{contig}~{post_contig_string}.bcf
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
        File output_bcf = "~{output_basename}.~{contig}~{post_contig_string}.bcf"
        File output_bcf_index = "~{output_basename}.~{contig}~{post_contig_string}.bcf.csi"
    }
}
