version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow BeagleImputationValidation {
    input {
        File ref_panel_vcf
        File ref_panel_vcf_index
        File truth_vcf
        File truth_vcf_index
        File test_vcf
        File test_vcf_index
        File? sample_names_to_select
        String output_basename
    }

    call RunBeagleGtStats {
        input:
            ref_panel_vcf = ref_panel_vcf,
            ref_panel_vcf_index = ref_panel_vcf_index,
            output_basename = output_basename
    }

    call SelectVariantType as SelectSnps {
        input:
            truth_vcf = truth_vcf,
            test_vcf = test_vcf,
            select_type_string = "SNP",
            sample_names_to_select = sample_names_to_select
    }

    call SelectVariantType as SelectIndels {
        input:
            truth_vcf = truth_vcf,
            test_vcf = test_vcf,
            select_type_string = "INDEL",
            sample_names_to_select = sample_names_to_select
    }

    call RunBeagleImputedR2 as RunBeagleImputedR2Snps {
        input:
            gt_stats = RunBeagleGtStats.gt_stats_output,
            truth_vcf = SelectSnps.truth_output_vcf,
            truth_vcf_index = SelectSnps.truth_output_vcf_index,
            test_vcf = SelectSnps.test_output_vcf,
            test_vcf_index = SelectSnps.test_output_vcf_index,
            output_basename = output_basename + ".SNPs"
    }

    call RunBeagleImputedR2 as RunBeagleImputedR2Indels {
        input:
            gt_stats = RunBeagleGtStats.gt_stats_output,
            truth_vcf = SelectIndels.truth_output_vcf,
            truth_vcf_index = SelectIndels.truth_output_vcf_index,
            test_vcf = SelectIndels.test_output_vcf,
            test_vcf_index = SelectIndels.test_output_vcf_index,
            output_basename = output_basename + ".INDELS"
    }

    output {
        File gt_stats_output = RunBeagleGtStats.gt_stats_output
        File imputed_r2_output_snps = RunBeagleImputedR2Snps.imputed_r2_output
        File imputed_r2_output_indels = RunBeagleImputedR2Indels.imputed_r2_output
    }
}

task RunBeagleGtStats {
    input {
        File ref_panel_vcf
        File ref_panel_vcf_index
        String output_basename

        Int disk_size_gb = ceil(size(ref_panel_vcf, "GiB"))  + 20
        Int cpu = 16
        Int memory_mb = 96000
    }

    Int command_mem = memory_mb - 6000
    Int max_heap = memory_mb - 4000

    command {
        set -e -o pipefail

        gunzip -c ~{ref_panel_vcf} | java -Xms~{command_mem}m -Xmx~{max_heap}m -jar /beagle_jars/gt-stats.jar > ~{output_basename}_gt_stats.tsv

    }

    output {
        File gt_stats_output = "~{output_basename}_gt_stats.tsv"
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/jsotobroad/beagle_validation:1.0.1"
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }
}

task RunBeagleImputedR2 {
    input {
        File gt_stats
        File truth_vcf
        File truth_vcf_index
        File test_vcf
        File test_vcf_index
        String output_basename

        Int disk_size_gb = ceil(size(truth_vcf, "GiB")) + ceil(size(test_vcf, "GiB")) + 20
        Int cpu = 16
        Int memory_mb = 100000
    }

    Int command_mem = memory_mb - 12000
    Int max_heap = memory_mb - 10000

    command {
        set -e -o pipefail

        java -Xms~{command_mem}m -Xmx~{max_heap}m -jar /beagle_jars/imputed-r2.jar ~{gt_stats} \
        ~{truth_vcf} \
        ~{test_vcf} > ~{output_basename}.imputed_stats.tsv

    }

    output {
        File imputed_r2_output = "~{output_basename}.imputed_stats.tsv"
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/jsotobroad/beagle_validation:1.0.1"
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }
}

task SelectVariantType {
    input {
        File truth_vcf
        File test_vcf
        String select_type_string = "SNP"
        File? sample_names_to_select

        Int disk_size_gb = ceil(2 * (size(truth_vcf, "GiB") + size(test_vcf, "GiB"))) + 10
        Int cpu = 1
        Int memory_mb = 6000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000
    String truth_basename = basename(truth_vcf, ".vcf.gz")
    String test_basename = basename(test_vcf, ".vcf.gz")

    command {
        set -euo pipefail

        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        SelectVariants \
        -select-type ~{select_type_string} \
        ~{"--sample-name " + sample_names_to_select} \
        --preserve-alleles \
        -V ~{truth_vcf} \
        -O ~{truth_basename}_~{select_type_string}.vcf.gz

        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        SelectVariants \
        -select-type ~{select_type_string} \
        ~{"--sample-name " + sample_names_to_select} \
        --preserve-alleles \
        -V ~{test_vcf} \
        -O ~{test_basename}_~{select_type_string}.vcf.gz
    }
    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }
    parameter_meta {
        truth_vcf: {
                      description: "truth vcf",
                      localization_optional: true
                  }
        test_vcf: {
                      description: "test vcf",
                      localization_optional: true
                  }
    }
    output {
        File truth_output_vcf = "~{truth_basename}_~{select_type_string}.vcf.gz"
        File truth_output_vcf_index = "~{truth_basename}_~{select_type_string}.vcf.gz.tbi"
        File test_output_vcf = "~{test_basename}_~{select_type_string}.vcf.gz"
        File test_output_vcf_index = "~{test_basename}_~{select_type_string}.vcf.gz.tbi"
    }
}
