version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow ExtractSamplesAndFilter {
    input {
        File input_bcf
        File input_bcf_index
        File sample_list
        String contig
        String output_basename
        String? post_contig_string
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

    output {
        File output_bcf = ExtractAndFilter.output_bcf
        File output_bcf_index = ExtractAndFilter.output_bcf_index
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
        # TODO change this when official docker image is ready
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:sps_sv_docker_images"
    }

    command <<<
        set -euo pipefail

        bcftools view \
            -S ~{sample_list} \                        # subset samples first
            -O b -o sample_subset.bcf

        bcftools view \
            -i 'GT[*]="alt" && INFO/AF >= 0.05' \      # keep alt sites (i.e. remove hom ref sites) and filter for AF
            sample_subset.bcf \
            -O b \                                     # export compressed bcf
            -o ~{output_basename}.~{contig}~{post_contig_string}.bcf
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
