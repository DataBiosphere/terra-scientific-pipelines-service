version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow MakeSitesOnly {
    input {
        File input_bcf
        File input_bcf_index
        String contig
        String output_basename
        String? post_contig_string
    }

    call DropGenotypes {
        input:
            input_bcf = input_bcf,
            input_bcf_index = input_bcf_index,
            contig = contig,
            output_basename = output_basename,
            post_contig_string = post_contig_string
    }

    output {
        File sites_only_bcf = DropGenotypes.output_bcf
        File sites_only_bcf_index = DropGenotypes.output_bcf_index
    }
}

task DropGenotypes {
    input {
        File input_bcf
        File input_bcf_index
        String contig
        String output_basename
        String post_contig_string = ""

        Int disk_size_gb = ceil(2 * (size(input_bcf, "GiB") + size(input_bcf_index, "GiB"))) + 20
        Int cpu = 1
        Int memory_mb = 12000
        # TODO change this when official docker image is ready
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:sps_sv_docker_images"
    }

    String output_basename_full =  "~{output_basename}.~{contig}~{post_contig_string}"

    command <<<
        set -euo pipefail

        # drop genotype (sample) columns, keeping sites only
        bcftools view -G ~{input_bcf} -Ob -o ~{output_basename_full}.bcf

        bcftools index ~{output_basename_full}.bcf
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
        File output_bcf = "~{output_basename_full}.bcf"
        File output_bcf_index = "~{output_basename_full}.bcf.csi"
    }
}
