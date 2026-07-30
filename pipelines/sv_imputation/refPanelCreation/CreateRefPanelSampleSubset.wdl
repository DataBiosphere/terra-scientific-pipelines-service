version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow CreateRefPanelSampleSubset {
    input {
        # from source ref panel
        File input_panel_bubble_bcf
        File input_panel_bubble_bcf_index
        File input_panel_id_split_vcf
        File input_panel_id_split_vcf_index

        # for generated panel
        File sample_list

        String contig
        String output_basename
        String? copy_to_cloud_dest
    }

    # extract samples from parent ref panel, drop hom ref sites
    call ExtractAndFilter {
        input:
            input_bcf = input_panel_bubble_bcf,
            input_bcf_index = input_panel_bubble_bcf_index,
            sample_list = sample_list,
            output_basename = "~{output_basename}.~{contig}"
    }

    # split to biallelic bcf
    call SplitMultiallellics {
        input:
            bcf_input = ExtractAndFilter.output_bcf,
            bcf_input_index = ExtractAndFilter.output_bcf_index,
            output_basename = "~{output_basename}.~{contig}.split"
     }

    # generate sites-only split bcf
    call MakeSitesOnly {
        input:
            input_bcf = SplitMultiallellics.output_bcf,
            input_bcf_index = SplitMultiallellics.output_bcf_index,
            output_basename = "~{output_basename}.~{contig}.split.sites"
    }

    # generate panel ids vcf
    call ExtractIdsAndFilter {
        input:
            split_sites_only_bcf = MakeSitesOnly.output_bcf,
            split_sites_only_bcf_index = MakeSitesOnly.output_bcf_index,
            input_panel_id_split_vcf = input_panel_id_split_vcf,
            input_panel_id_split_vcf_index = input_panel_id_split_vcf_index,
            output_basename = "~{output_basename}.~{contig}.id.split"
    }

    output {
        File bubble_split_bcf = SplitMultiallellics.output_bcf
        File bubble_split_bcf_index = SplitMultiallellics.output_bcf_index
        File bubble_split_sites_bcf = MakeSitesOnly.output_bcf
        File bubble_split_sites_bcf_index = MakeSitesOnly.output_bcf_index
        File panel_id_split_vcf = ExtractIdsAndFilter.output_panel_id_split_vcf
        File panel_id_split_vcf_index = ExtractIdsAndFilter.output_panel_id_split_vcf_index
    }
}

task ExtractAndFilter {
    input {
        File input_bcf
        File input_bcf_index
        File sample_list
        String output_basename

        Int disk_size_gb = ceil(3 * (size(input_bcf, "GiB") + size(input_bcf_index, "GiB"))) + 20
        Int cpu = 1
        Int memory_mb = 12000
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
    }

    command <<<
        set -euo pipefail

        # subset samples -> drop hom ref sites
        bcftools view -S ~{sample_list} -Ou ~{input_bcf} \
            | bcftools view -i 'GT[*]="alt"' - -Ob -o ~{output_basename}.bcf

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


task SplitMultiallellics {
    input {
        File bcf_input
        File bcf_input_index
        String output_basename

        Int disk_size_gb =  ceil(3 * (size(bcf_input, "GiB") + size(bcf_input_index, "GiB"))) + 20
        Int cpu = 1
        Int memory_mb = 24000
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
    }


    command <<<
        set -e -o pipefail

        # split multiallelics w/o normalization -> recalculate AC, AN -> drop AC=0
        bcftools norm -m -any -N ~{bcf_input} -Ou | \
            bcftools +fill-tags - -Ou -- -t AC,AN | \
            bcftools view --min-ac 1 -Ob -o "~{output_basename}.bcf"

        bcftools index "~{output_basename}.bcf"
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


task MakeSitesOnly {
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

task ExtractIdsAndFilter {
    input {
        File split_sites_only_bcf
        File split_sites_only_bcf_index
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

        # extract a unique list of all INFO/ID values; a single record's INFO/ID field may
        # contain multiple IDs separated by colons, so split those onto
        # their own lines before deduplicating. drop empty and missing (".") entries.
        bcftools query -f '%INFO/ID\n' ~{split_sites_only_bcf} \
            | tr ':' '\n\n' \
            | sed -e '/^$/d' -e '/^\.$/d' \
            | sort -u > ids_to_keep.list

        # keep only records whose INFO/ID is in the ids_to_keep list
        bcftools view \
            -i 'INFO/ID=@ids_to_keep.list' \
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
