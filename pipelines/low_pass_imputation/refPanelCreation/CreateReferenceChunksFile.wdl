version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow CreateReferenceChunksFile {
    input {
        Array[Array[File]] split_reference_chunks_files
        Array[String] contigs
        String output_basename

        # if you want to use the memory values (column 3 of output) from a pre-existing file, provide that file here
        File? source_for_memory_values_file
        Int? default_memory_value_gb
    }

    Array[String] split_reference_chunks_files_flat = flatten(split_reference_chunks_files)

    call CreateRefChunksFile {
        input:
            split_reference_chunks_files = split_reference_chunks_files_flat,
            contigs = contigs,
            output_basename = output_basename,
            source_for_memory_values_file = source_for_memory_values_file,
            default_memory_value_gb = default_memory_value_gb
    }

    output {
        File reference_chunks_file = CreateRefChunksFile.output_reference_chunks_file
    }
}

task CreateRefChunksFile {
    input {
        Array[String] split_reference_chunks_files
        Array[String] contigs
        String output_basename
        File? source_for_memory_values_file
        Int default_memory_value_gb = 16
    }

    command {
        set -e -o pipefail

        # the input file is just a list of the .bin file paths

        # the output file will have three columns: the first two are the chromosome name (format "chr1") and the path to the
        # .bin chunk file. The third is the memory value to use for that chunk when running GLIMPSE2.

        # if source_for_memory_values_file is not provided, all memory values will be set to the value of default_memory_value_gb.
        # if source_for_memory_values_file is provided, it should be a tab-delimited file with three columns, of which
        # the third is the memory values in GB for each chunk. The first column should exactly match the first column
        # of the output file.

        default_memory_value_gb="~{default_memory_value_gb}"

        printf "%s\n" ~{sep='\n' split_reference_chunks_files} > /tmp/reference_chunks_files.txt

        if [ -n "~{source_for_memory_values_file}" ]; then
            cat "~{source_for_memory_values_file}" | grep ~{sep='|' contigs} | > /tmp/memory_values_file.tsv
        fi

        cat <<'EOF' > script.py

import pandas as pd
import os
import sys

reference_chunks_files = pd.read_csv("/tmp/reference_chunks_files.txt", header=None, names=["chunk_file"])
if "~{source_for_memory_values_file}" != "":
    memory_values = pd.read_csv("/tmp/memory_values_file.tsv", sep="\t", header=None, names=["chrom", "chunk_file", "memory_gb"])




        EOF
        python3 script.py
    }

    output {
        File output_reference_chunks_file = "~{output_basename}.reference_chunks.txt"
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
    }
}