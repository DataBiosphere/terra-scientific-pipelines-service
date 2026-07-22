version 1.0

import "CopyToCloud.wdl" as CopyToCloudTask

# This script is under review. It is not actively tested or maintained at this time.
workflow SplitMultiallelicsBcf {
    input {
        File input_bcf
        File input_bcf_index
        File ref_dict
        String contig
        String output_basename
        Int num_base_chunk_size = 10000000
        String? copy_to_cloud_dest
    }

    String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"

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
                bcf = input_bcf,
                bcf_index = input_bcf_index,
                start = start_chunk_first,
                end = end_chunk_first,
                chrom = contig,
                basename = chunk_basename_first
        }

        call SeparateMultiallelics {
            input:
                bcf_input = GenerateChunk.output_bcf,
                bcf_input_index = GenerateChunk.output_bcf_index,
        }
    }

    call GatherBcfs {
        input:
            input_bcfs = SeparateMultiallelics.output_bcf,
            output_bcf_name = output_basename + ".multiallelic_split." + contig + ".bcf"
    }

    if (defined(copy_to_cloud_dest)) {
        call CopyToCloudTask.CopyToCloud as CopyToCloud {
            input:
                source_file = GatherBcfs.output_bcf,
                source_file_index = GatherBcfs.output_bcf_index,
                copy_to_cloud_dest = select_first([copy_to_cloud_dest])
        }
    }

    output {
        String multi_allelics_split_bcf = select_first([CopyToCloud.copied_file, GatherBcfs.output_bcf])
        String multi_allelics_split_bcf_index = select_first([CopyToCloud.copied_file_index, GatherBcfs.output_bcf_index])
    }
}

task SeparateMultiallelics {
    input {
        File bcf_input
        File bcf_input_index

        Int disk_size_gb =  ceil(3 * (size(bcf_input, "GiB") + size(bcf_input_index, "GiB"))) + 20
        Int cpu = 1
        Int memory_mb = 12000
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
    }

    String output_basename = basename(bcf_input, ".bcf")

    command <<<
        set -e -o pipefail

        bcftools norm -m - ~{bcf_input} -O b -o ~{output_basename}.multiallelic_split.bcf

        bcftools index ~{output_basename}.multiallelic_split.bcf
    >>>

    output {
        File output_bcf = "~{output_basename}.multiallelic_split.bcf"
        File output_bcf_index = "~{output_basename}.multiallelic_split.bcf.csi"
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

task CalculateChromosomeLength {
    input {
        File ref_dict
        String chrom

        String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
        Int memory_mb = 2000
        Int cpu = 1
        Int disk_size_gb = ceil(2 * size(ref_dict, "GiB")) + 5
    }

    command <<<
        set -e -o pipefail

        grep -P "SN:~{chrom}\t" ~{ref_dict} | sed 's/.*LN://' | sed 's/\t.*//'
    >>>

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
        File bcf
        File bcf_index

        Int disk_size_gb = ceil(3 * (size(bcf, "GiB") + size(bcf_index, "GiB"))) + 20
        Int cpu = 1
        Int memory_mb = 6000
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
    }

    command <<<
        set -euo pipefail

        # extract the region for this chunk, excluding filtered records
        bcftools view \
            -r ~{chrom}:~{start}-~{end} \
            -f PASS,. \
            ~{bcf} \
            -O b -o ~{basename}.bcf

        bcftools index ~{basename}.bcf
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
        File output_bcf = "~{basename}.bcf"
        File output_bcf_index = "~{basename}.bcf.csi"
    }
}

task GatherBcfs {
    input {
        Array[File] input_bcfs
        String output_bcf_name

        Int disk_size_gb = ceil(1.5 * size(input_bcfs, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 6000
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
    }

    command <<<
        set -euo pipefail

        # the chunks are non-overlapping and in genomic order, so they can be concatenated directly
        bcftools concat ~{sep=" " input_bcfs} -O b -o ~{output_bcf_name}

        bcftools index ~{output_bcf_name}
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
        File output_bcf = "~{output_bcf_name}"
        File output_bcf_index = "~{output_bcf_name}.csi"
    }
}
