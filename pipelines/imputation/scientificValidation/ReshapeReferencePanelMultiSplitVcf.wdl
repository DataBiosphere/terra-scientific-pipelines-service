version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow ReshapeReferencePanelSplitVcf {
    input {
        File ref_panel_vcf
        File ref_panel_vcf_index
        String output_base_name
        Boolean use_bcftools
        File ref_dict
        String contig
        Int num_base_chunk_size
        Int sample_chunk_size

        String ubuntu_docker = "ubuntu:20.04"
        String gatk_docker = "broadinstitute/gatk:4.6.0.0"
    }

    call ChunkSampleNames {
        input:
            vcf = ref_panel_vcf,
            sample_chunk_size = sample_chunk_size
    }

    Float sample_chunk_size_float = sample_chunk_size
    Int num_sample_chunks = ceil(ChunkSampleNames.sample_count / sample_chunk_size_float)

    call CalculateChromosomeLength {
        input:
            ref_dict = ref_dict,
            chrom = contig,
            ubuntu_docker = ubuntu_docker
    }

    Float num_base_chunk_float = num_base_chunk_size
    Int num_base_chunks = ceil(CalculateChromosomeLength.chrom_length / num_base_chunk_float)

    #scatter (i in range(num_base_chunks)) {
    scatter (i in range(2)) {
        Int start_chunk_first = (i * num_base_chunk_size) + 1
        Int end_chunk_first = if (CalculateChromosomeLength.chrom_length < ((i + 1) * num_base_chunk_size)) then CalculateChromosomeLength.chrom_length else ((i + 1) * num_base_chunk_size)
        String chunk_basename_first = "generate_first_chunk_" + i

        call GenerateChunk as GenerateChunkFirst {
            input:
                vcf = ref_panel_vcf,
                vcf_index = ref_panel_vcf_index,
                start = start_chunk_first,
                end = end_chunk_first,
                chrom = contig,
                basename = chunk_basename_first,
                gatk_docker = gatk_docker
        }

        scatter (j in range(num_sample_chunks)) {
            Int start_sample = (j * sample_chunk_size) + 10
            Int end_sample = if (ChunkSampleNames.sample_count <= ((j + 1) * sample_chunk_size)) then ChunkSampleNames.sample_count + 9 else ((j + 1) * sample_chunk_size ) + 9

            call SelectSamplesWithCut {
                input:
                    vcf = GenerateChunkFirst.output_vcf,
                    cut_start_field = start_sample,
                    cut_end_field = end_sample,
                    chunk_index = j
            }
        }
    }
    Array[Array[File]] vcfs_to_be_gathered_first = transpose(SelectSamplesWithCut.output_vcf)

    scatter (i in range(length(vcfs_to_be_gathered_first))) {
        call GatherVcfs as GatherVcfsFirst {
            input:
                input_vcfs = vcfs_to_be_gathered_first[i],
                output_vcf_name = "gather_vcf_first_chunk_" + i
        }
    }

    # run reshape here on gathervcfsfirst output

    scatter(i in range(length(GatherVcfsFirst.output_vcf))) {
        #scatter (j in range(num_base_chunks)) {
        scatter (j in range(2)) {
            Int start_chunk_second = (j * num_base_chunk_size) + 1
            Int end_chunk_second = if (CalculateChromosomeLength.chrom_length < ((j + 1) * num_base_chunk_size)) then CalculateChromosomeLength.chrom_length else ((j + 1) * num_base_chunk_size)
            String chunk_second_basename = "generate_second_chunk_" + j

            call GenerateChunk as GenerateChunkSecond {
                input:
                    vcf = GatherVcfsFirst.output_vcf[i],
                    vcf_index = GatherVcfsFirst.output_vcf_index[i],
                    start = start_chunk_second,
                    end = end_chunk_second,
                    chrom = contig,
                    basename = chunk_second_basename,
                    gatk_docker = gatk_docker
            }
        }
    }

    Array[Array[File]] vcfs_to_be_gathered_second = transpose(GenerateChunkSecond.output_vcf)

    scatter (i in range(length(vcfs_to_be_gathered_second))) {
        call MergeVcfsWithCutPaste {
            input:
                vcfs = vcfs_to_be_gathered_second[i],
                basename = "merge_chunk_" + i
        }
    }

    call GatherVcfs as GatherVcfsSecond {
        input:
            input_vcfs = MergeVcfsWithCutPaste.output_vcf,
            output_vcf_name = output_base_name
    }

    output {
        File output_vcf = GatherVcfsSecond.output_vcf
        File output_vcf_index = GatherVcfsSecond.output_vcf_index
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

task SelectSamplesWithCut {
    input {
        File vcf

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

        mkfifo fifo_bgzip
        mkfifo fifo_cut

        bcftools view -h --no-version ~{vcf} | awk '!/^#CHROM/' > header.vcf
        n_lines=$(wc -l header.vcf | cut -d' ' -f1)

        cat header.vcf
        echo $n_lines

        bgzip -d ~{vcf} -o fifo_bgzip &
        tail +$((n_lines)) fifo_bgzip | cut -f 1-9,~{cut_start_field}-~{cut_end_field} > fifo_cut &

        cat header.vcf fifo_cut | bgzip -o ~{basename(vcf)}.sample_chunk_~{chunk_index}.vcf.gz

        bctools index -t ~{basename(vcf)}.sample_chunk_~{chunk_index}.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/ckachulis/bcftools_bgzip"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: memory_mb + " MiB"
        cpu: cpu
    }

    output {
        File output_vcf = "~{basename(vcf)}.chunk_~{chunk_index}.vcf.gz"
        File output_vcf_index = "~{basename(vcf)}.chunk_~{chunk_index}.vcf.gz.tbi"
    }

}

task MergeVcfsWithCutPaste {
    input {
        Array[File] vcfs
        String basename

        Int disk_size_gb = ceil(2.2 * size(vcfs, "GiB") + 10)
        Int mem_gb = 10
        Int cpu = 10
        Int preemptible = 0
    }

    command <<<
        set -euo pipefail

        vcfs=(~{sep=" " vcfs})

        mkfifo fifo_0
        mkfifo fifo_to_paste_0

        i=1

        fifos_to_paste=()
        bcftools view -h --no-version ${vcfs[0]} | awk '!/^#CHROM/' > header.vcf
        n_lines=$(wc -l header.vcf | cut -d' ' -f1)

        cat header.vcf
        echo $n_lines

        bgzip -d ${vcfs[0]} -o fifo_0 &

        tail +$((n_lines+1)) fifo_0 > fifo_to_paste_0 &

        for vcf in "${vcfs[@]:1}"; do
        fifo_name="fifo_$i"
        mkfifo "$fifo_name"

        fifo_name_to_paste="fifo_to_paste_$i"
        mkfifo "$fifo_name_to_paste"
        fifos_to_paste+=("$fifo_name_to_paste")
        n_lines=$(bcftools view -h --no-version $vcf | awk '!/^#CHROM/' | wc -l | cut -d' ' -f1)

        bgzip -d ${vcf} -o "$fifo_name" &
        tail +$((n_lines+1)) "$fifo_name" | cut -f 10- > "$fifo_name_to_paste" &

        ((i++))
        done

        mkfifo fifo_to_cat

        paste fifo_to_paste_0 "${fifos_to_paste[@]}" > fifo_to_cat &

        cat header.vcf fifo_to_cat | bgzip -o ~{basename}.merged.vcf.gz

        bctools index -t ~{basename}.merged.vcf.gz
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
        File output_vcf_index = "~{basename}.merged.vcf.gz.tbi"
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

task CalculateChromosomeLength {
    input {
        File ref_dict
        String chrom

        String ubuntu_docker = "ubuntu:20.04"
        Int memory_mb = 2000
        Int cpu = 1
        Int disk_size_gb = ceil(2*size(ref_dict, "GiB")) + 5
    }

    command {
        set -e -o pipefail

        grep -P "SN:~{chrom}\t" ~{ref_dict} | sed 's/.*LN://' | sed 's/\t.*//'
    }
    runtime {
        docker: ubuntu_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 3
    }
    output {
        Int chrom_length = read_int(stdout())
    }
}

task GenerateChunk {
    input {
        Int start
        Int end
        String chrom
        String basename
        File vcf
        File vcf_index

        Int disk_size_gb = ceil(2*size(vcf, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 8000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    command {
        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        SelectVariants \
        -V ~{vcf} \
        -L ~{chrom}:~{start}-~{end} \
        -O ~{basename}.vcf.gz \
        --exclude-filtered true
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
        File output_vcf = "~{basename}.vcf.gz"
        File output_vcf_index = "~{basename}.vcf.gz.tbi"
    }
}

task GatherVcfs {

    input {
        Array[File] input_vcfs
        String output_vcf_name
        Int disk_size_gb = ceil(1.5 * size(input_vcfs, "GiB")) + 10
        Int machine_mem_mb = 6000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    }

    parameter_meta {
        input_vcfs: {
                        localization_optional: true
                    }
    }

    command <<<
        set -euo pipefail

        # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
        # This argument disables expensive checks that the file headers contain the same set of
        # genotyped samples and that files are in order by position of first record.
        gatk --java-options "-Xms4000m -Xmx4500m" \
        GatherVcfsCloud \
        --ignore-safety-checks \
        --gather-type BLOCK \
        --input ~{sep=" --input " input_vcfs} \
        --output ~{output_vcf_name}

        tabix ~{output_vcf_name}
    >>>

    runtime {
        memory: "~{machine_mem_mb} MiB"
        cpu: "1"
        bootDiskSizeGb: 15
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
        docker: gatk_docker
    }

    output {
        File output_vcf = "~{output_vcf_name}"
        File output_vcf_index = "~{output_vcf_name}.tbi"
    }
}
