version 1.0


workflow CreatePanelAuxiliaryFiles {
    input {
        Array[String] chromosomes

        # per chromosome, in same order, used to create pop_glimpse2_panel_resources_json
        Array[File] panel_bubble_split_sites_only_vcf_array
        Array[File] panel_bubble_split_sites_only_vcf_idx_array
        Array[File] panel_id_split_vcf_gz_array
        Array[File] panel_id_split_vcf_gz_tbi_array

        # these files are gathered together to make preprocess_panel_bubble_split_sites_only_vcf input
        Array[File] reduced_panel_bubble_split_simple_sites_bcf
        Array[File] reduced_panel_bubble_split_simple_sites_bcf_index

        File ref_dict
        File reference_fasta_fai
        String output_prefix

        Int interval_size = 25000000
    }

    scatter (i in range(length(chromosomes))) {
        String chromosome = chromosomes[i]

        call CalculateChromosomeLength {
            input:
                ref_dict = ref_dict,
                chrom = chromosome
        }

        call GenerateChromosomeIntervals {
            input:
                chromosome   = chromosome,
                chrom_length = CalculateChromosomeLength.chrom_length,
                interval_size = interval_size
        }

        PopAndMarginalizePanelResourcesChromosome pop_panel_resources_chromosome = object {
            panel_bubble_split_sites_only_vcf: panel_bubble_split_sites_only_vcf_array[i],
            panel_bubble_split_sites_only_vcf_idx: panel_bubble_split_sites_only_vcf_idx_array[i],
            panel_id_split_vcf_gz: panel_id_split_vcf_gz_array[i],
            panel_id_split_vcf_gz_tbi: panel_id_split_vcf_gz_tbi_array[i],
            pop_regions: GenerateChromosomeIntervals.intervals
        }
        Pair[String, PopAndMarginalizePanelResourcesChromosome] pop_panel_resources_chromosome_pair = (chromosome, pop_panel_resources_chromosome)
    }

    call CoercePairsToMap {
        input:
            pair_array = pop_panel_resources_chromosome_pair,
            output_prefix = output_prefix
    }

    call GatherAndSortBcfs {
        input:
            bcf_files = reduced_panel_bubble_split_simple_sites_bcf,
            bcf_index_files = reduced_panel_bubble_split_simple_sites_bcf_index,
            reference_fasta_fai = reference_fasta_fai,
            output_basename = output_prefix + ".reduced_panel_bubble_split_simple_sites"
    }

    output {
        File pop_panel_resources_json = CoercePairsToMap.out_map_json
        File gathered_reduced_panel_bubble_split_simple_sites_bcf = GatherAndSortBcfs.sorted_bcf
        File gathered_reduced_panel_bubble_split_simple_sites_bcf_index = GatherAndSortBcfs.sorted_bcf_index
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Boolean? use_ssd
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

struct PopAndMarginalizePanelResourcesChromosome {
    String panel_bubble_split_sites_only_vcf
    String panel_bubble_split_sites_only_vcf_idx
    String panel_id_split_vcf_gz
    String panel_id_split_vcf_gz_tbi
    Array[String] pop_regions
}

task CoercePairsToMap {
    input {
        Array[Pair[String, PopAndMarginalizePanelResourcesChromosome]] pair_array
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        python3 <<CODE
        import json

        # write_json is safe here because it executes inside the task container
        with open("~{write_json(pair_array)}", "r") as f:
            pairs = json.load(f)

        out_map = {item["left"]: item["right"] for item in pairs}

        with open("~{output_prefix}.json", "w") as f:
            json.dump(out_map, f, indent=2)
        CODE
    >>>

    output {
        File out_map_json = "~{output_prefix}.json"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            10,
        boot_disk_gb:       10,
        use_ssd:            false,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-pypy:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CalculateChromosomeLength {
    input {
        File ref_dict
        String chrom

        String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
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

task GenerateChromosomeIntervals {
    input {
        String chromosome
        Int chrom_length
        Int interval_size

        Int memory_mb = 2000
        Int cpu = 1
        Int disk_size_gb = 10
        String docker = "us.gcr.io/broad-dsp-gcr-public/base/python:3.12-debian"
    }

    command <<<
        python3 <<CODE
        chrom = "~{chromosome}"
        length = ~{chrom_length}
        step = ~{interval_size}

        intervals = []
        start = 1
        while start <= length:
            end = min(start + step - 1, length)
            intervals.append(f"{chrom}:{start}-{end}")
            start = end + 1

        with open("intervals.txt", "w") as f:
            f.write("\n".join(intervals) + "\n")
        CODE
    >>>

    output {
        Array[String] intervals = read_lines("intervals.txt")
    }

    runtime {
        docker: docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
    }
}

task GatherAndSortBcfs {
    input {
        Array[File] bcf_files
        Array[File] bcf_index_files
        File reference_fasta_fai
        String output_basename

        Int cpu = 4
        Int memory_mb = 16000
        Int disk_size_gb = ceil(4 * size(bcf_files, "GiB")) + 20
        String docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
    }

    command <<<
        set -e -o pipefail

        # Sort the input BCF file list by FAI chromosome order before concatenating.
        # Each input BCF already covers exactly one chromosome and is internally sorted,
        # so concatenating in FAI order (chr1, chr2, chr3 ...) is sufficient —
        # no global bcftools sort step is needed or used.
        python3 <<CODE
import subprocess, sys

fai_order = {}
with open("~{reference_fasta_fai}") as f:
    for i, line in enumerate(f):
        fai_order[line.split('\t')[0].strip()] = i

bcf_paths = [l.strip() for l in open("~{write_lines(bcf_files)}") if l.strip()]

def get_chrom(path):
    r = subprocess.run(
        ['bcftools', 'query', '-f', '%CHROM\n', path],
        capture_output=True, text=True, check=True
    )
    lines = [l.strip() for l in r.stdout.split('\n') if l.strip()]
    if not lines:
        sys.exit(f"ERROR: no records found in {path}")
    return lines[0]

sorted_paths = sorted(bcf_paths, key=lambda p: fai_order.get(get_chrom(p), 999999))

with open('sorted_filelist.txt', 'w') as f:
    f.write('\n'.join(sorted_paths) + '\n')
CODE

        # Concatenate in FAI-sorted order
        bcftools concat \
            --file-list sorted_filelist.txt \
            --allow-overlaps \
            -Ob -o concat.bcf
        bcftools index concat.bcf

        # Rebuild ##contig lines in FAI order in the header.
        # bcftools concat inherits ##contig lines from the first input file's header,
        # which may still be in alphabetical order. ##contig lines must come before
        # the #CHROM line or tools will report "CHROM is not defined".
        bcftools view -h concat.bcf | grep -v "^##contig" | grep -v "^#CHROM" > new_header.txt
        awk '{print "##contig=<ID="$1",length="$2">"}' ~{reference_fasta_fai} >> new_header.txt
        bcftools view -h concat.bcf | grep "^#CHROM" >> new_header.txt

        bcftools reheader -h new_header.txt concat.bcf -o ~{output_basename}.sorted.bcf
        bcftools index ~{output_basename}.sorted.bcf
    >>>

    output {
        File sorted_bcf       = "~{output_basename}.sorted.bcf"
        File sorted_bcf_index = "~{output_basename}.sorted.bcf.csi"
    }

    runtime {
        docker:  docker
        disks:   "local-disk ${disk_size_gb} HDD"
        memory:  "${memory_mb} MiB"
        cpu:     cpu
        noAddress: true
    }
}

