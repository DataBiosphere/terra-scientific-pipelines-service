version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow CreateImputationRefPanelBeagle {
    input {
        File ref_vcf
        File ref_vcf_index
        String chromosome
        Boolean create_brefs = true
        Boolean create_interval_lists = true
        Boolean create_bed_files = true
        File ref_dict
        # this is used to chunk up the input vcfs to create interval lists from them in a timely manner
        Int chunkLength = 2500000
        String output_basename
    }

    Float chunkLengthFloat = chunkLength

    String custom_basename_with_chr = output_basename + "." + chromosome

    call CalculateChromosomeLength {
        input:
            ref_dict = ref_dict,
            chrom = chromosome
    }

    Int num_chunks = ceil(CalculateChromosomeLength.chrom_length / chunkLengthFloat)

    if (create_interval_lists || create_bed_files) {
        scatter (i in range(num_chunks)) {
            String custom_basename_with_chr_and_chunk = output_basename + "." + chromosome + ".chunk_" + i

            Int start = (i * chunkLength) + 1
            Int end = if (CalculateChromosomeLength.chrom_length < ((i + 1) * chunkLength)) then CalculateChromosomeLength.chrom_length else ((i + 1) * chunkLength)

            call CreateRefPanelIntervalLists {
                input:
                    ref_panel_vcf = ref_vcf,
                    ref_panel_vcf_index = ref_vcf_index,
                    output_basename = custom_basename_with_chr_and_chunk,
                    chrom = chromosome,
                    start = start,
                    end = end
            }
        }

        call GatherIntervalLists as GatherChunkedIntervalLists {
            input:
                basename = custom_basename_with_chr,
                interval_lists = CreateRefPanelIntervalLists.interval_list
        }
    }

    if (create_bed_files) {
        call CreateRefPanelBedFiles {
            input:
                ref_panel_interval_list = select_first([GatherChunkedIntervalLists.interval_list]),
                basename = custom_basename_with_chr,
        }
    }

    if (create_brefs) {
        call BuildBref3 {
            input:
                vcf = ref_vcf,
                basename = custom_basename_with_chr
        }
    }

    output {
        Array[File?] interval_lists = select_first([GatherChunkedIntervalLists.interval_list, []])
        Array[File?] bed_files = select_first([CreateRefPanelBedFiles.bed_file, []])
        Array[File?] brefs = select_first([BuildBref3.bref3, []])
    }
}

task CreateRefPanelIntervalLists {
    input {
        File ref_panel_vcf
        File ref_panel_vcf_index

        String chrom
        Int start
        Int end

        String? output_basename

        Int disk_size_gb = ceil(size(ref_panel_vcf, "GiB") / 2) # not sure how big the disk size needs to be since we aren't downloading the entire VCF here
        Int cpu = 1
        Int memory_mb = 6000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    }

    Int command_mem = memory_mb - 2500
    Int max_heap = memory_mb - 2000

    String name_from_file = basename(ref_panel_vcf, ".vcf.gz")
    String basename = select_first([output_basename, name_from_file])

    command {
        set -e -o pipefail

        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        SelectVariants \
        -V ~{ref_panel_vcf} \
        -O ~{basename}.vcf.gz \
        -L ~{chrom}:~{start}-~{end}

        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        VcfToIntervalList \
        -I ~{basename}.vcf.gz \
        -O ~{basename}.interval_list \
        --MAX_RECORDS_IN_RAM 100000
    }

    output {
        File interval_list = "~{basename}.interval_list"
    }

    parameter_meta {
        ref_panel_vcf: {
                           description: "vcf",
                           localization_optional: true
                       }
        ref_panel_vcf_index: {
                                 description: "vcf index",
                                 localization_optional: true
                             }
    }

    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
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
    }
    output {
        Int chrom_length = read_int(stdout())
    }
}

task GatherIntervalLists {
    input {
        Array[File] interval_lists

        String basename

        Int disk_size_gb = 50
        Int cpu = 1
        Int memory_mb = 12000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    }

    Int command_mem = memory_mb - 2500
    Int max_heap = memory_mb - 2000

    command {
        set -e -o pipefail

        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        IntervalListTools \
        ACTION=CONCAT \
        SORT=true \
        UNIQUE=true \
        I=~{sep=' I= ' interval_lists} \
        O=~{basename}.interval_list
    }

    output {
        File interval_list = "~{basename}.interval_list"
    }

    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }
}

task CreateRefPanelBedFiles {
    input {
        File ref_panel_interval_list
        String basename
        Int disk_size_gb = ceil(2*size(ref_panel_interval_list, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 12000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    }

    Int command_mem = memory_mb - 2000
    Int max_heap = memory_mb - 1500


    command {
        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        IntervalListToBed \
        -I ~{ref_panel_interval_list} \
        -O ~{basename}.bed
    }

    output {
        File bed_file = "~{basename}.bed"
    }

    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }
}

task BuildBref3 {
    input {
        File vcf
        String basename
        # this can be smaller if you have a small input vcf but we do have a large gap between command mem
        # and machine mem
        Int memory_mb = 50000
        # bref outputs are very small, 50 gb overhead is plenty
        Int disk_size =  ceil(size(vcf, "GiB")) + 50
    }

    Int command_mem = memory_mb - 10000
    Int max_heap = memory_mb - 5000

    command <<<
        java -Xms~{command_mem}m -Xmx~{max_heap}m -jar /usr/gitc/bref3.jar ~{vcf} > ~{basename}.bref3
    >>>

    runtime {
        docker: "jsotobroad/beagle_jar:bit_shift_2"
        memory: "${memory_mb} MiB"
        cpu: 4
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File bref3 = "~{basename}.bref3"
    }
}
