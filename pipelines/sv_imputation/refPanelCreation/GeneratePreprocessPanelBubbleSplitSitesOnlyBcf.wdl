version 1.0
    
workflow GeneratePreprocessPanelBubbleSplitSitesOnlyBcf {
    input {
        File reduced_panel_bubble_vcf
        File reduced_panel_bubble_vcf_index
        String output_basename
        String contig
    }

    call SelectSimpleSites {
        input:
            reduced_panel_bubble_vcf = reduced_panel_bubble_vcf,
            reduced_panel_bubble_vcf_index = reduced_panel_bubble_vcf_index,
            output_basename = output_basename + "." + contig + ".simple",
    }

    call CreateBcfIndex {
        input:
            bcf_input = SelectSimpleSites.output_bcf
    }

    output {
        File output_bcf = CreateBcfIndex.output_bcf
        File output_bcf_index = CreateBcfIndex.output_bcf_index
    }
}

task CreateBcfIndex {
    input {
        File bcf_input

        Int disk_size_gb = ceil(1.2 * size(bcf_input, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 6000
    }

    String bcf_basename = basename(bcf_input)

    command {
        set -e -o pipefail

        ln -sf ~{bcf_input} ~{bcf_basename}

        bcftools index ~{bcf_basename}
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-pypy:v1"
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        noAddress: true
    }

    output {
        File output_bcf = "~{bcf_basename}"
        File output_bcf_index = "~{bcf_basename}.csi"
    }
}

task SelectSimpleSites {
    input {
        File reduced_panel_bubble_vcf
        File reduced_panel_bubble_vcf_index
        String output_basename

        Int disk_size_gb = ceil(4 * size(reduced_panel_bubble_vcf, "GiB")) + 30
        Int memory_mb = 24000
    }

    command {
        set -e -o pipefail

        echo "doing norm"
        bcftools norm -m -any ~{reduced_panel_bubble_vcf} -Ou | bcftools view -G -Ob -o multialleic_split.bcf
        echo "doing query"
        bcftools query -f '%INFO/ID\n' multialleic_split.bcf | grep -v ":" > simple.ids.list
        echo "doing view"
        bcftools view -i "INFO/ID=@simple.ids.list" -Ob -o ~{output_basename}.bcf
    }

    output {
        File output_bcf = "~{output_basename}.bcf"
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-pypy:v1"
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        preemptible: 0
        cpu: 1
        noAddress: true
    }
}
