version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow ReshapeReferencePanel {
    input {
        File ref_panel_vcf
        File genetic_map
        String chromosome
        String output_basename
        Int reshape_threads
    }

    call ReshapeReferencePanel {
        input:
            ref_panel_vcf = ref_panel_vcf,
            genetic_map = genetic_map,
            chrom = chromosome,
            output_basename = output_basename,
            num_threads = reshape_threads

    }

    output {
        File recombined_reference_panel = ReshapeReferencePanel.output_vcf
    }
}

task ReshapeReferencePanel {
    input {
        File ref_panel_vcf
        File genetic_map
        String chrom
        String output_basename
        Int num_threads

        Int disk_size_gb = ceil(3*size(ref_panel_vcf, "GiB")) + 20
        Int memory_mb = 6000
    }

    command {
        set -e -o pipefail

        haploshuffling_static --vcf ~{ref_panel_vcf} \
        --region ~{chrom} \
        --output ~{output_basename}.~{chrom}.reshaped.vcf.gz \
        --map ~{genetic_map} \
        --seed 12345 \
        --gen 8
        --threads ~{num_threads}

    }

    output {
        File output_vcf = "~{output_basename}.~{chrom}.reshaped.vcf.gz"
    }

    runtime {
        docker: "theocavinato/reshape"
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: num_threads
    }
}
