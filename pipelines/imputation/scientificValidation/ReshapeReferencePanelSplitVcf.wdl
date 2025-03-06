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

    call ConvertVcfToBcf {
        input:
            vcf = ref_panel_vcf
    }

#    call ChunkSampleNames {
#        input:
#            vcf = ref_panel_vcf_header,
#            sample_chunk_size = sample_chunk_size
#    }
#
#    Float sample_chunk_size_float = sample_chunk_size
#    Int num_chunks = ceil(ChunkSampleNames.sample_count / sample_chunk_size_float)
#
#    #scatter (i in range(num_chunks)) {
#    scatter (i in range(1)) {
#        Int start = (i * sample_chunk_size) + 10
#        Int end = if (ChunkSampleNames.sample_count <= ((i + 1) * sample_chunk_size)) then ChunkSampleNames.sample_count + 9 else ((i + 1) * sample_chunk_size ) + 9
##        if (localize_vcfs && use_gatk) {
##            call SelectSamplesFromVcfWithGatkLocalize {
##                input:
##                    vcf = ref_panel_vcf,
##                    vcf_index = ref_panel_vcf_index,
##                    monitoring_script = monitoring_script,
##                    sample_names = ChunkSampleNames.sample_names[i],
##                    chunk_index = i
##            }
##        }
##
##        if (!localize_vcfs && use_gatk) {
##            call SelectSamplesFromVcfWithGatkStream {
##                input:
##                    vcf = ref_panel_vcf,
##                    vcf_index = ref_panel_vcf_index,
##                    monitoring_script = monitoring_script,
##                    sample_names = ChunkSampleNames.sample_names[i],
##                    chunk_index = i
##            }
##        }
##
##        if (!use_gatk) {
##            call SelectSamplesFromVcfWithBcftools {
##                input:
##                    vcf = ref_panel_vcf,
##                    vcf_index = ref_panel_vcf_index,
##                    sample_names = ChunkSampleNames.sample_names[i],
##                    chunk_index = i
##            }
##        }
##
##        File select_output = select_first([SelectSamplesFromVcfWithGatkLocalize.output_vcf, SelectSamplesFromVcfWithGatkStream.output_vcf, SelectSamplesWithCut.output_vcf])
#        call SelectSamplesWithCut {
#            input:
#                vcf = ref_panel_vcf,
#                monitoring_script = monitoring_script,
#                cut_start_field = start,
#                cut_end_field = end,
#                chunk_index = i
#        }
#    }
#
#    call MergeVcfsWithCutPaste {
#        input:
#            vcfs = SelectSamplesWithCut.output_vcf,
#            monitoring_script = monitoring_script,
#            basename = output_base_name
#    }
#
#    call CreateVcfIndex {
#        input:
#            vcf_input = MergeVcfsWithCutPaste.output_vcf
#    }
#
#    output {
#        File reshaped_reference_panel = CreateVcfIndex.output_vcf
#        File reshaped_reference_panel_index = CreateVcfIndex.output_vcf_index
#    }
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

        # store number of samples in a file to be used in downstream scattering/tasks
        wc -l sample_names.txt | cut -d' ' -f1 > sample_count.txt
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
        Int sample_count = read_int("sample_count.txt")
    }
}

task ConvertVcfToBcf {
    input {
        File vcf

        String basename = basename(vcf, '.vcf.gz')

        Int disk_size_gb = ceil(size(vcf, "GiB")) + 10
        String bcftools_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/denovo:2025-02-11-v1.0.2-hotfix-22bf77e0"
        Int cpu = 1
        Int memory_mb = 6000
    }
    command {
        set -e -o pipefail

        bcftools view -Ob ~{vcf} > ~{basename}.bcf
        bcftools index ~{basename}.bcf
    }

    runtime {
        docker: bcftools_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }

    output {
        File sample_names = "~{basename}.bcf"
        File sample_names = "~{basename}.bcf.csi"
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

task SelectSamplesWithCut {
    input {
        File vcf
        File monitoring_script

        Int cut_start_field
        Int cut_end_field
        Int chunk_index

        Int disk_size_gb = ceil(1.5*size(vcf, "GiB")) + 10
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int cpu = 2
        Int memory_mb = 6000
    }

    command <<<
        set -euo pipefail

        bash ~{monitoring_script} &

        mkfifo fifo_bgzip
        mkfifo fifo_cut

        bcftools view -h --no-version ~{vcf} | awk '!/^#CHROM/' > header.vcf
        n_lines=$(wc -l header.vcf | cut -d' ' -f1)

        cat header.vcf
        echo $n_lines

        bgzip -d ~{vcf} -o fifo_bgzip &
        tail +$((n_lines)) fifo_bgzip | cut -f 1-9,~{cut_start_field}-~{cut_end_field} > fifo_cut &

        cat header.vcf fifo_cut | bgzip -o ~{basename(vcf)}.chunk_~{chunk_index}.vcf.gz

        bcftools view -h ~{basename(vcf)}.chunk_~{chunk_index}.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/ckachulis/bcftools_bgzip"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: memory_mb + " MiB"
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
        File output_vcf = "~{vcf_basename}"
        File output_vcf_index = "~{vcf_basename}.tbi"
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

task MergeVcfsWithCutPaste {
    input {
        Array[File] vcfs
        File monitoring_script
        String basename

        Int disk_size_gb = ceil(2.2 * size(vcfs, "GiB") + 10)
        Int mem_gb = 16
        Int cpu = 2
        Int preemptible = 0
    }

    command <<<
        set -euo pipefail

        bash ~{monitoring_script} &

        vcfs=(~{sep=" " vcfs})

        mkfifo fifo_0
        mkfifo fifo_to_paste_0

        i=1

        fifos_to_paste=()
        bcftools view -h --no-version ${vcfs[0]} | awk '!/^#CHROM/' > header.vcf
        n_lines=$(wc -l header.vcf | cut -d' ' -f1)

        bgzip -d ${vcfs[0]} -o fifo_0 &

        tail +$((n_lines+1)) fifo_0 > fifo_to_paste_0 &

        for vcf in "${vcfs[@]:1}"; do
            fifo_name="fifo_$i"
            mkfifo "$fifo_name"

            fifo_name_to_paste="fifo_to_paste_$i"
            mkfifo "$fifo_name_to_paste"
            fifos_to_paste+=("$fifo_name_to_paste")
            n_lines=$(bcftools view -h --no-version $vcf | wc -l | cut -d' ' -f1)

            bgzip -d ${vcf} -o "$fifo_name" &
            tail +$((n_lines)) "$fifo_name" | cut -f 10- > "$fifo_name_to_paste" &

            ((i++))
        done

        mkfifo fifo_to_cat

        paste fifo_to_paste_0 "${fifos_to_paste[@]}" > fifo_to_cat &

        cat header.vcf fifo_to_cat | bgzip -o ~{basename}.merged.vcf.gz

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/ckachulis/bcftools_bgzip"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File output_vcf = "~{basename}.merged.vcf.gz"
    }

}
