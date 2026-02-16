version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow LeftAlignVcf {
    input {
        File input_vcf
        File input_vcf_index
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        String contig
        Int max_indel_length = 200
        Int num_base_chunk_size = 10000000
    }

    String basename = basename(input_vcf, ".vcf.gz")

    String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"

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

        call BcftoolsNorm {
            input:
                vcf = GenerateChunk.output_vcf,
                vcf_index = GenerateChunk.output_vcf_index,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                basename = basename + ".bcftools." + contig
        }
    }

    call GatherVcfs as GatherBcftoolsVcfs {
        input:
            input_vcfs = BcftoolsNorm.output_vcf,
            output_vcf_name = basename + "." + contig + ".left_aligned.vcf.gz"
    }

    output {
        File left_aligned_vcf = GatherBcftoolsVcfs.output_vcf
        File left_aligned_vcf_index = GatherBcftoolsVcfs.output_vcf_index
    }
}

task BcftoolsNorm {
    input {
        String basename
        File vcf
        File vcf_index
        File ref_fasta
        File ref_fasta_index
        File ref_dict

        Int disk_size_gb = ceil(1.5 * size(vcf, "GiB") + size(ref_fasta, "GiB") + size(ref_dict, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 6000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    command {
        set -euo pipefail

        bcftools norm -f ~{ref_fasta} ~{vcf} -o ~{basename}.left_aligned.vcf.gz -Oz

        bcftools index -t ~{basename}.left_aligned.vcf.gz
    }
    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }

    output {
        File output_vcf = "~{basename}.left_aligned.vcf.gz"
        File output_vcf_index = "~{basename}.left_aligned.vcf.gz.tbi"
    }
}

task GatherVcfs {
    input {
        Array[File] input_vcfs
        String output_vcf_name
        Int disk_size_gb = ceil(1.2 * size(input_vcfs, "GiB")) + 10
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

task CreateVcfIndex {
    input {
        File vcf_input

        Int disk_size_gb = ceil(1.2 * size(vcf_input, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 6000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
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
        Int memory_mb = 6000
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
