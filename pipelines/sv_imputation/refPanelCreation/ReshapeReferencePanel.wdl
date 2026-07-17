version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow ReshapeReferencePanel {
    input {
        File ref_panel_bcf
        String output_basename
        File genetic_map
        String contig
        Int reshape_threads

        String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.0.0"
    }

    call ReshapeReferencePanel {
        input:
            ref_panel_bcf = ref_panel_bcf,
            genetic_map = genetic_map,
            chrom = contig,
            output_basename = output_basename,
            num_threads = reshape_threads
    }

        call CreateBcfIndex as CreateBcfIndexReshapeReferencePanel {
            input:
                bcf_input = ReshapeReferencePanel.output_bcf
        }

    output {
        File output_vcf = CreateBcfIndexReshapeReferencePanel.output_bcf
        File output_vcf_index = CreateBcfIndexReshapeReferencePanel.output_bcf_index
    }
}

task CreateBcfIndex {
    input {
        File bcf_input

        Int disk_size_gb = ceil(1.2 * size(bcf_input, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 6000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    }

    String bcf_basename = basename(bcf_input)

    command {
        set -e -o pipefail

        ln -sf ~{bcf_input} ~{bcf_basename}

        bcftools index ~{bcf_basename}
    }

    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }

    output {
        File output_bcf = "~{bcf_basename}"
        File output_bcf_index = "~{bcf_basename}.csi"
    }
}

task ReshapeReferencePanel {
    input {
        File ref_panel_bcf
        File genetic_map
        String chrom
        String output_basename
        Int num_threads

        Int disk_size_gb = ceil(3 * size(ref_panel_bcf, "GiB")) + 20
        Int memory_mb = 6000
    }

    command {
        set -e -o pipefail

        haploshuffling_static --vcf ~{ref_panel_bcf} \
        --region ~{chrom} \
        --output ~{output_basename}.~{chrom}.reshaped.bcf \
        --map ~{genetic_map} \
        --seed 12345 \
        --gen 8 \
        --threads ~{num_threads}

    }

    output {
        File output_bcf = "~{output_basename}.~{chrom}.reshaped.bcf"
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/jsotobroad/theocavinato/reshape"
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        preemptible: 0
        cpu: num_threads
    }
}
