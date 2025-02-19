version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow ReshapeReferencePanelSplitVcf {
    input {
        File ref_panel_vcf
        Int sample_chunk_size = 50000
    }

    call ChunkSampleNames {
        input:
            vcf = ref_panel_vcf,
            sample_chunk_size = sample_chunk_size
    }

    scatter (i in range(length(ChunkSampleNames.sample_name_args))) {
        call SelectSamplesFromVcfWithGatk {
            input:
                vcf = ref_panel_vcf,
                sample_name_args = ChunkSampleNames.sample_name_args[i],
                chunk_index = i
        }
    }

    output {
        Array[File] recombined_reference_panel = SelectSamplesFromVcfWithGatk.output_vcf
    }
}

task ChunkSampleNames {
    input {
        File vcf
        Int sample_chunk_size

        Int disk_size_gb = ceil(size(vcf, "GiB")) + 10
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int cpu = 1
        Int memory_mb = 4000
    }
    command <<<
        set -e -o pipefail

        # query the sample names from the VCF and chunk by sample_chunk_size
        bcftools query -l ~{vcf} > sample_names.txt
        split -l ~{sample_chunk_size} sample_names.txt sample_chunk_args

        # add the prefix "--sample-name" to each line in file to use as a GATK arguments file
        for f in sample_chunk_args*; do sed -i 's/^/--sample-name /' $f; done

    >>>
    output {
        Array[File] sample_name_args = glob("sample_chunk_args*")
    }
    runtime {
        docker: bcftools_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }
}


task SelectSamplesFromVcfWithGatk {
    input {
        File vcf
        File sample_name_args
        Int chunk_index

        Int disk_size_gb = 800 # streaming so unsure how big this vcf is going to be
        Int cpu = 1
        Int memory_mb = 6000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    }

    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    String basename = basename(vcf, '.vcf.gz')

    command {
        set -e -o pipefail

        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        SelectVariants \
        -V ~{vcf} \
        --arguments_file ~{sample_name_args} \
        -O ~{basename}.chunk_~{chunk_index}.vcf.gz
    }

    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }

    parameter_meta {
        vcf: {
                 description: "vcf",
                 localization_optional: true
             }
    }
    output {
        File output_vcf = "~{basename}.chunk_~{chunk_index}.vcf.gz"
    }
}
