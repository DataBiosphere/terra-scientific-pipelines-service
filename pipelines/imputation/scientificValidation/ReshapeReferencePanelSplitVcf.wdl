version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow ReshapeReferencePanelSplitVcf {
    input {
        File ref_panel_vcf
        File ref_panel_vcf_index
        File ref_panel_vcf_header # this is possibly created during ref panel creation
        File monitoring_script
        String output_base_name
        Boolean localize_vcfs
        Boolean use_gatk
        Int sample_chunk_size
    }

    call ChunkSampleNames {
        input:
            vcf = ref_panel_vcf_header,
            sample_chunk_size = sample_chunk_size
    }

    scatter (i in range(length(ChunkSampleNames.sample_names))) {
        if (localize_vcfs && use_gatk) {
            call SelectSamplesFromVcfWithGatkLocalize {
                input:
                    vcf = ref_panel_vcf,
                    vcf_index = ref_panel_vcf_index,
                    monitoring_script = monitoring_script,
                    sample_names = ChunkSampleNames.sample_names[i],
                    chunk_index = i
            }
        }

        if (!localize_vcfs && use_gatk) {
            call SelectSamplesFromVcfWithGatkStream {
                input:
                    vcf = ref_panel_vcf,
                    vcf_index = ref_panel_vcf_index,
                    monitoring_script = monitoring_script,
                    sample_names = ChunkSampleNames.sample_names[i],
                    chunk_index = i
            }
        }

        if (!use_gatk) {
            call SelectSamplesFromVcfWithBcftools {
                input:
                    vcf = ref_panel_vcf,
                    vcf_index = ref_panel_vcf_index,
                    sample_names = ChunkSampleNames.sample_names[i],
                    chunk_index = i
            }
        }

        File select_output = select_first([SelectSamplesFromVcfWithGatkLocalize.output_vcf, SelectSamplesFromVcfWithGatkStream.output_vcf, SelectSamplesFromVcfWithBcftools.output_vcf])

        call CreateVcfIndex {
            input:
                vcf_input = select_output
        }
    }

    call MergeVcfsBcfTools {
        input:
            input_vcfs = CreateVcfIndex.vcf,
            input_vcf_indices = CreateVcfIndex.vcf_index,
            output_vcf_basename = output_base_name
    }

    output {
        File recombined_reference_panel = MergeVcfsBcfTools.output_vcf
    }
}

task ChunkSampleNames {
    input {
        File vcf
        Int sample_chunk_size

        Int disk_size_gb = ceil(size(vcf, "GiB")) + 10
        String bcftools_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/denovo:2025-02-11-v1.0.2-hotfix-22bf77e0"
        Int cpu = 1
        Int memory_mb = 4000
    }
    command {
        set -e -o pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        # query the sample names from the VCF and chunk by sample_chunk_size
        bcftools head ~{vcf} > vcf_header.txt
        bcftools query -l vcf_header.txt > sample_names.txt
        split -l ~{sample_chunk_size} sample_names.txt sample_chunks
    }

    runtime {
        docker: bcftools_docker
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
        Array[File] sample_names = glob("sample_chunks*")
    }
}


task SelectSamplesFromVcfWithGatkStream {
    input {
        File vcf
        File vcf_index
        File sample_names
        File monitoring_script
        Int chunk_index

        Int disk_size_gb = ceil(0.5*size(vcf, "GiB")) + 10
        Int cpu = 2
        Int memory_mb = 16000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    }

    Int command_mem = memory_mb - 6000
    Int max_heap = memory_mb - 4000

    String basename = basename(vcf, '.vcf.gz')

    command {
        set -e -o pipefail

        bash ~{monitoring_script} &

        # add the prefix "--sample-name" to each line in file to use as a GATK arguments file
        cp ~{sample_names} sample_names_gatk_args.txt
        sed -i 's/^/--sample-name /' sample_names_gatk_args.txt

        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        SelectVariants \
        -V ~{vcf} \
        --arguments_file sample_names_gatk_args.txt \
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
        vcf_index: {
                 description: "vcf index",
                 localization_optional: true
             }
    }

    output {
        File output_vcf = "~{basename}.chunk_~{chunk_index}.vcf.gz"
    }
}

task SelectSamplesFromVcfWithGatkLocalize {
    input {
        File vcf
        File vcf_index
        File sample_names
        File monitoring_script
        Int chunk_index

        Int disk_size_gb = ceil(1.5*size(vcf, "GiB")) + 10
        Int cpu = 2
        Int memory_mb = 16000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    }

    Int command_mem = memory_mb - 6000
    Int max_heap = memory_mb - 4000

    String basename = basename(vcf, '.vcf.gz')

    command {
        set -e -o pipefail

        bash ~{monitoring_script} &

        # add the prefix "--sample-name" to each line in file to use as a GATK arguments file
        cp ~{sample_names} sample_names_gatk_args.txt
        sed -i 's/^/--sample-name /' sample_names_gatk_args.txt

        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        SelectVariants \
        -V ~{vcf} \
        --arguments_file sample_names_gatk_args.txt \
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
                 localization_optional: false
             }
        vcf_index: {
                       description: "vcf index",
                       localization_optional: false
                   }
    }

    output {
        File output_vcf = "~{basename}.chunk_~{chunk_index}.vcf.gz"
    }
}

task SelectSamplesFromVcfWithBcftools {
    input {
        File vcf
        File vcf_index
        File sample_names
        Int chunk_index

        Int disk_size_gb = ceil(1.5*size(vcf, "GiB")) + 10
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int cpu = 1
        Int memory_mb = 6000
    }

    command {
        set -e -o pipefail

        # query the sample names from the VCF and chunk by sample_chunk_size
        bcftools view -S ~{sample_names} -O z -o ~{basename(vcf)}.chunk_~{chunk_index}.vcf.gz ~{vcf}
    }

    runtime {
        docker: bcftools_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }

    output {
        File output_vcf = "~{basename(vcf)}.chunk_~{chunk_index}.vcf.gz"
    }
}

task CreateVcfIndex {
    input {
        File vcf_input

        Int disk_size_gb = ceil(1.2*size(vcf_input, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 6000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    }

    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    String vcf_basename = basename(vcf_input)

    command {
        set -e -o pipefail

        ln -sf ~{vcf_input} ~{vcf_basename}

        bcftools index -t ~{vcf_basename}
    }

    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }

    output {
        File vcf = "~{vcf_basename}"
        File vcf_index = "~{vcf_basename}.tbi"
    }
}

task MergeVcfsBcfTools {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        String output_vcf_basename

        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int memory_mb = 10000
        Int cpu = 3
        Int disk_size_gb = 3 * ceil(size(input_vcfs, "GiB") + size(input_vcf_indices, "GiB")) + 20
    }

    command {
        set -e -o pipefail

        bcftools merge --threads ~{cpu} ~{sep=' ' input_vcfs} -O z -o ~{output_vcf_basename}.vcf.gz
    }

    runtime {
        docker: bcftools_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }

    output {
        File output_vcf = "~{output_vcf_basename}.vcf.gz"
    }
}
