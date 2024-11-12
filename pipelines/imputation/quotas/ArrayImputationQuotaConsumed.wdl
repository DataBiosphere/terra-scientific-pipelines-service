version 1.0

workflow QuotaConsumed {
    input {
        Int chunkLength = 25000000
        Int chunkOverlaps = 5000000

        File multi_sample_vcf

        File ref_dict
        Array[String] contigs
        String reference_panel_path_prefix
        String genetic_maps_path
        String output_basename
        Boolean split_output_to_single_sample = false

        # file extensions used to find reference panel files
        String interval_list_suffix = ".interval_list"
        String bref3_suffix = ".bref3"
    }

    call CountSamples {
        input:
            vcf = multi_sample_vcf
    }

    output {
        Int quota_consumed = CountSamples.nSamples
    }
}

task CountSamples {
    input {
        File vcf

        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int cpu = 1
        Int memory_mb = 3000
        Int disk_size_gb = ceil(size(vcf, "GiB")) + 10
    }

    command <<<
        bcftools query -l ~{vcf} | wc -l
    >>>
    runtime {
        docker: bcftools_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }
    output {
        Int nSamples = read_int(stdout())
    }
}
