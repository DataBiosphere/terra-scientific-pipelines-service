version 1.0

workflow Glimpse2SplitReference {
    input {
        # There are two variables required to set define the contigs. The first one, contig_regions, defines the coordinates for the separate chunks,
        # usually this would be an array of ["chr1", "chr2", ..., "chr22", "chrX:1-2781479", "chrX:2781480-155701382", "chrX:155701383-156030895"].
        # Notice the split of chrX into PAR1, non-PAR, PAR2. When constructing the reference panel chunks, however, we need to determine the appropriate
        # reference panel file names. We could do that by parsing the substring before the colon (which WDL doesn't support), or we can just pass
        # another array contig_names_in_reference_panel = ["chr1", "chr2", ..., "chr22", "chrX", "chrX", "chrX"]. Note that contig_regions and
        # contig_names_in_reference_panel must have the same length.
        String contig_name
        File contig_reference_chunks

        # Example path for chr1: gs://bucket/to/panel/reference_panel_chr1_merged.bcf(.csi)
        # reference_panel_prefix = "gs://bucket/to/panel/reference_panel_"
        # reference_panel_suffix = "_merged.bcf"
        # reference_panel_index_suffix = ".csi"
        # String reference_panel_prefix
        # String reference_panel_suffix
        # String reference_panel_index_suffix
        String reference_filename
        String reference_filename_index

        # Same format as reference_panel_pre-/suffix
        String genetic_map_path_prefix
        String genetic_map_path_suffix

        File ref_dict

        Int seed = 42
        Boolean keep_monomorphic_ref_sites = true

        Boolean add_allele_info = true

        Int? ac_cutoff

        Int generate_chunk_default_memory_mb = 6000
        Int fix_annotations_default_memory_gb = 6
        Int glimpse_default_memory_gb = 16
        String? generate_chunk_memory_override
        String? fix_annotations_memory_override
        String? glimpse_split_reference_memory_override

        # New docker same as old but with bcftools v1.21 instead of v1.16
        String docker = "us.gcr.io/broad-dsde-methods/updated_glimpse_docker:v1.0"
    }

    String generate_chunk_memory_override_defined = select_first([generate_chunk_memory_override, "{'empty': 'empty}"])
    String fix_annotations_memory_override_defined = select_first([fix_annotations_memory_override, "{'empty': 'empty}"])
    String glimpse_split_reference_memory_override_defined = select_first([glimpse_split_reference_memory_override, "{'empty': 'empty}"])


    File json_file = write_json(generate_chunk_memory_override_defined)

    # String reference_filename = reference_panel_prefix + contig_name + reference_panel_suffix
    # String reference_filename_index = reference_filename + reference_panel_index_suffix
    String genetic_map_filename = genetic_map_path_prefix + contig_name + genetic_map_path_suffix

    # Shard the VCF file into chunks and process for GLIMPSE
    Array[String] contig_reference_chunks_lines = read_lines(contig_reference_chunks)

    call ConvertJsonStringToMap as GenerateChunkConvertString {
        input:
            json_string = generate_chunk_memory_override_defined
    }

    call BuildMemoryMap as GenerateChunkMemoryMap {
        input:
            memory_override_map = read_json(json_file),
            default_memory_gb = generate_chunk_default_memory_mb,
            num_shards = length(contig_reference_chunks_lines)
    }

    call ConvertJsonStringToMap as FixAnnotationsConvertString {
        input:
            json_string = fix_annotations_memory_override_defined
    }

    call BuildMemoryMap as FixAnnotationsMemoryMap {
        input:
            memory_override_map = FixAnnotationsConvertString.output_map,
            default_memory_gb = fix_annotations_default_memory_gb,
            num_shards = length(contig_reference_chunks_lines)
    }

    call ConvertJsonStringToMap as GlimpseConvertString {
        input:
            json_string = glimpse_split_reference_memory_override_defined
    }

    call BuildMemoryMap as GlimpseMemoryMap {
        input:
            memory_override_map = GlimpseConvertString.output_map,
            default_memory_gb = glimpse_default_memory_gb,
            num_shards = length(contig_reference_chunks_lines)
    }

    call CalculateChromosomeLength {
        input:
            ref_dict = ref_dict,
            chrom = contig_name,
            ubuntu_docker = docker
    }

    scatter (i in range(length(contig_reference_chunks_lines))) {
        String interval = contig_reference_chunks_lines[i]
        String generate_chunk_memory_value = GenerateChunkMemoryMap.memory_values[i]
        String fix_annotations_memory_value = FixAnnotationsMemoryMap.memory_values[i]
        String glimpse_memory_value = GlimpseMemoryMap.memory_values[i]

        call GenerateChunk {
            input:
                vcf = reference_filename,
                vcf_index = reference_filename_index,
                interval = interval,
                contig_length = CalculateChromosomeLength.chrom_length,
                memory_mb = generate_chunk_memory_value
        }

        if (add_allele_info) {
            call FixAnnotations {
                input:
                    vcf = GenerateChunk.output_vcf,
                    vcf_index = GenerateChunk.output_vcf_index,
                    interval = interval,
                    mem_gb = fix_annotations_memory_value
            }
        }

        File sharded_vcf = select_first([FixAnnotations.vcf_chunk, GenerateChunk.output_vcf])
        File sharded_vcf_index = select_first([FixAnnotations.vcf_chunk_index, GenerateChunk.output_vcf_index])

        # Apply allele count cutoff if specified
        if (defined(ac_cutoff)) {
            call ApplyACCutoff {
                input:
                    vcf = sharded_vcf,
                    ac_cutoff = select_first([ac_cutoff]),
            }
        }

        call GlimpseSplitReferenceTask {
            input:
                reference_panel = select_first([ApplyACCutoff.filtered_vcf, sharded_vcf]),
                reference_panel_index = select_first([ApplyACCutoff.filtered_vcf_index, sharded_vcf_index]),
                contig = contig_name,
                interval = interval,
                genetic_map = genetic_map_filename,
                seed = seed,
                keep_monomorphic_ref_sites = keep_monomorphic_ref_sites,
                docker = docker,
                mem_gb = glimpse_memory_value
        }
    }


    output {
        Array[File] reference_chunks = flatten(GlimpseSplitReferenceTask.split_reference_chunks)
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

task ConvertJsonStringToMap {
    input {
        String json_string

        String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
        Int memory_mb = 2000
        Int cpu = 1
        Int disk_size_gb = 10
    }

    command {
        set -e -o pipefail

        echo ~{json_string} >
    }
    runtime {
        docker: ubuntu_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 3
    }
    output {
        Map[String, String] output_map = read_json(stdout())
    }
}

task BuildMemoryMap {
    input {
        Map[String, String] memory_override_map = {"empty": "empty"} # Default to empty map if not provided
        Int default_memory_gb
        Int num_shards

        # Runtime parameters
        Int disk_size = 10
        Int preemptible = 3
    }

    command <<<
        python <<CODE
        import json

        # Get inputs
        memory_override_file = "~{write_json(select_first([memory_override_map]))}"
        with open(memory_override_file, "r") as f:
            memory_override = json.load(f)
        default_memory = ~{default_memory_gb}
        num_shards = ~{num_shards}

        # Create memory map
        memory_map = {}
        for i in range(num_shards):
            if str(i) in memory_override:
                memory_map[str(i)] = memory_override[str(i)]
            else:
                memory_map[str(i)] = str(default_memory)

        # Write the memory map to output file
        # with open("memory_map.json", "w") as f:
        #     json.dump(memory_map, f)

        # Write array of memory values to output file
        memory_values = list(memory_map.values())
        with open('memory_values.txt', 'w') as f:
            f.writelines('\n'.join(memory_values))

        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        memory: "2 GiB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible
    }

    output {
        Array[String] memory_values = read_lines("memory_values.txt")
    }
}

task GenerateChunk {
    input {
        String interval
        Int contig_length
        File vcf
        File vcf_index

        Int disk_size_gb = ceil(0.25 * size(vcf, "GiB")) + 30
        Int cpu = 1
        Int memory_mb = 6000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    command <<<
        set -euo pipefail

        INTERVAL=$(echo "~{interval}" | awk '{ print $3 }')

        # the last interval given by glimpse is not the end of the contig
        # which caused gatk to fail, this check makes sure the max value
        # for `end` is at most the length of the contig

        IFS=':|-' read -r chrom start end <<< $INTERVAL
        end=$(( end < ~{contig_length} ? end : ~{contig_length} ))

        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        SelectVariants \
        -V ~{vcf} \
        -L $chrom:$start-$end \
        -O chunked.vcf.gz \
        --exclude-filtered true
    >>>
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
        File output_vcf = "chunked.vcf.gz"
        File output_vcf_index = "chunked.vcf.gz.tbi"
    }
}

task FixAnnotations {
    input {
        File vcf
        File vcf_index

        String interval

        Int mem_gb = 6
        Int disk_size = ceil(2.5 * size(vcf, "GiB") + 50)
    }

    command <<<
        set -xueo pipefail

        INTERVAL=$(echo "~{interval}" | awk '{ print $3 }')

        # Use bcftools to recalculate AC, AN, and AF annotations
        bcftools view -r $INTERVAL ~{vcf} --threads $(nproc) -Ou | bcftools norm -Ou -m -any - | bcftools +fill-tags - -o "chunk.vcf.gz" -Wtbi -- -t AC,AN,AF
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/updated_glimpse_docker:v1.0"
        memory: mem_gb + " GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

    output {
        File vcf_chunk = "chunk.vcf.gz"
        File vcf_chunk_index = "chunk.vcf.gz.tbi"
    }
}

task ApplyACCutoff {
    input {
        File vcf
        Int ac_cutoff

        Int mem_gb = 6
        Int disk_size_gb = ceil(2.5 * size(vcf, "GiB") + 50)
    }

    command <<<
        set -xueo pipefail

        bcftools view -i "AC>=~{ac_cutoff}" ~{vcf} --threads $(nproc) -o "filtered.vcf.gz" -Wtbi
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/updated_glimpse_docker:v1.0"
        memory: mem_gb + " GB"
        cpu: 1
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }

    output {
        File filtered_vcf = "filtered.vcf.gz"
        File filtered_vcf_index = "filtered.vcf.gz.tbi"
    }
}

task GlimpseSplitReferenceTask {
    input {
        String contig
        File reference_panel
        File reference_panel_index
        File genetic_map
        String interval

        Int? seed
        Boolean keep_monomorphic_ref_sites

        String mem_gb = 16
        Int cpu = 4
        Int disk_size_gb = ceil(2.2 * size(reference_panel, "GiB") + size(genetic_map, "GiB") + 50)
        String docker
    }

    String keep_monomorphic_ref_sites_string = if keep_monomorphic_ref_sites then "--keep-monomorphic-ref-sites" else ""

    command <<<
        set -xeuo pipefail

        NPROC=$(nproc)
        echo "nproc reported ${NPROC} CPUs, using that number as the threads argument for GLIMPSE."

        # Print contig index to variable
        CONTIGINDEX="~{contig}"

        # Make chunk index from interval
        IN_INTERVAL=$(echo "~{interval}" | awk '{ print $3 }')
        OUT_INTERVAL=$(echo "~{interval}" | awk '{ print $4 }')
        CHUNKINDEX=$(echo "${IN_INTERVAL}" | tr ":" "-")

        mkdir -p reference_output_dir

        /bin/GLIMPSE2_split_reference \
            --threads ${NPROC} \
            --reference ~{reference_panel} \
            --map ~{genetic_map} \
            --input-region ${IN_INTERVAL} \
            --output-region ${OUT_INTERVAL} \
            --output reference_output_dir/reference_panel_contigindex_${CONTIGINDEX}_chunkindex_${CHUNKINDEX} \
            ~{keep_monomorphic_ref_sites_string} \
            ~{"--seed " + seed}
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: 0
    }

    output {
        Array[File] split_reference_chunks = glob("reference_output_dir/*.bin")
    }
}
