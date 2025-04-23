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
        File gt_stats_file
        String output_basename
    }

#    call RunBeagleGtStats {
#        input:
#            ref_panel_vcf = ref_panel_vcf,
#            ref_panel_vcf_index = ref_panel_vcf_index,
#            output_basename = output_basename
#    }

    call RunBeagleImputedR2 {
        input:
#            gt_stats = RunBeagleGtStats.gt_stats_output,
            gt_stats = gt_stats_file,
            truth_vcf = truth_vcf,
            truth_vcf_index = truth_vcf_index,
            test_vcf = test_vcf,
            test_vcf_index = test_vcf_index,
            output_basename = output_basename
    }

    output {
#        File gt_stats_output = RunBeagleGtStats.gt_stats_output
        File imputed_r2_output = RunBeagleImputedR2.imputed_r2_output
    }
}

task RunBeagleGtStats {
    input {
        File ref_panel_vcf
        File ref_panel_vcf_index
        String output_basename

        Int disk_size_gb = ceil(size(ref_panel_vcf, "GiB"))  + 20
        Int cpu = 8
        Int memory_mb = 63000
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
        docker: "us.gcr.io/broad-dsde-methods/jsotobroad/beagle_validation:1.0.0"
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
        Int cpu = 1
        Int memory_mb = 6000
    }

    Int command_mem = memory_mb - 3000
    Int max_heap = memory_mb - 2000

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
        docker: "us.gcr.io/broad-dsde-methods/jsotobroad/beagle_validation:1.0.0"
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }
}
