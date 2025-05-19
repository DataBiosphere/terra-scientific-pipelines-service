version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow SplitMultiallelics {
    input {
        File input_vcf
        File input_vcf_index
        File ref_dict
        String contig
        Int num_base_chunk_size = 10000000
    }

    String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.0.0"

    call CalculateChromosomeLength {
        input:
            ref_dict = ref_dict,
            chrom = contig,
            ubuntu_docker = ubuntu_docker
    }


    Float num_base_chunk_float = num_base_chunk_size
    Int num_base_chunks = ceil(CalculateChromosomeLength.chrom_length / num_base_chunk_float)

    scatter (i in range(num_base_chunks)) {
        Int start_chunk_first = (i * num_base_chunk_size) + 1
        Int end_chunk_first = if (CalculateChromosomeLength.chrom_length < ((i + 1) * num_base_chunk_size)) then CalculateChromosomeLength.chrom_length else ((i + 1) * num_base_chunk_size)
        String chunk_basename_first = "generate_first_chunk_" + i

        call GenerateChunk {
            input:
                vcf = input_vcf,
                vcf_index = input_vcf_index,
                start = start_chunk_first,
                end = end_chunk_first,
                chrom = contig,
                basename = chunk_basename_first,
                gatk_docker = gatk_docker
        }

        call SeparateMultiallelics {
            input:
                vcf_input = GenerateChunk.output_vcf,
                vcf_input_index = GenerateChunk.output_vcf_index,
        }
    }

    call GatherVcfs {
        input:
            input_vcfs = SeparateMultiallelics.output_vcf,
            output_vcf_name = basename(input_vcf, ".vcf.gz") + ".multiallelic_split.vcf.gz",
    }

    output {
        File multi_allelics_split_vcf = GatherVcfs.output_vcf
        File multi_allelics_split_vcf_index = GatherVcfs.output_vcf_index
    }
}

task SeparateMultiallelics {
    input {
        File vcf_input
        File vcf_input_index

        Int disk_size_gb =  ceil(2.5 * (size(vcf_input, "GiB") + size(vcf_input_index, "GiB"))) + 10
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int cpu = 1
        Int memory_mb = 6000
    }

    String output_basename = basename(vcf_input, ".vcf.gz")
    command {
        set -e -o pipefail

        bcftools norm -m - ~{vcf_input} -Oz -o ~{output_basename}.multiallelic_split.vcf.gz
    }
    output {
        File output_vcf = "~{output_basename}.multiallelic_split.vcf.gz"
    }
    runtime {
        docker: bcftools_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 0
        noAddress: true
    }
}

task CreateVcfIndex {
    input {
        File vcf_input

        Int disk_size_gb = ceil(1.2 * size(vcf_input, "GiB")) + 10
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

        String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
        Int memory_mb = 2000
        Int cpu = 1
        Int disk_size_gb = ceil(2 * size(ref_dict, "GiB")) + 5
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

        Int disk_size_gb = ceil(2 * size(vcf, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 8000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    command {
        set -euo pipefail

        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        SelectVariants \
        -V ~{vcf} \
        -L ~{chrom}:~{start}-~{end} \
        -O ~{basename}.vcf.gz \
        --exclude-filtered true \
        -select 'POS >= ~{start}'
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

    Int command_mem = machine_mem_mb - 1500
    Int max_heap = machine_mem_mb - 1000

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
        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
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
