version 1.0

workflow SubsetVcfByContigsAndSamples {
    input {
        File input_vcf
        File input_vcf_index
        File ref_fasta
        File ref_fasta_index
        File ref_fasta_dict
        Array[String] contigs
        File sample_names_file
        String output_basename

        Int preemptible = 0
    }

    scatter (contig in contigs) {
        call SubsetVcfBySampleNameAndContig {
            input:
                input_vcf = input_vcf,
                input_vcf_index = input_vcf_index,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_fasta_dict = ref_fasta_dict,
                contig = contig,
                sample_names_file = sample_names_file,
                output_basename = output_basename + "." + contig,
                preemptible = preemptible
        }
    }

    call ConcatenateVcfs {
        input:
            input_vcfs = SubsetVcfBySampleNameAndContig.output_vcf,
            input_vcf_indexes = SubsetVcfBySampleNameAndContig.output_vcf_index,
            output_basename = output_basename,
            preemptible = preemptible
    }

    output {
        File subset_vcf = ConcatenateVcfs.output_vcf
        File subset_vcf_index = ConcatenateVcfs.output_vcf_index
    }
}

# Subset the VCF down to the requested samples and one chromosome at a time.
task SubsetVcfBySampleNameAndContig {
    input {
        File input_vcf
        File input_vcf_index
        File ref_fasta
        File ref_fasta_index
        File ref_fasta_dict
        String contig
        File sample_names_file
        String output_basename

        Int disk_size_gb = ceil(3 * size(input_vcf, "GiB") + 20)
        Int cpu = 1
        Int memory_mb = 6000
        Int preemptible = 0
    }

    String output_vcf_name = output_basename + "." + contig + ".vcf.gz"

    command <<<
        set -euo pipefail

        gatk --java-options "-Xms~{memory_mb - 1500}m -Xmx~{memory_mb - 1000}m" \
            SelectVariants \
            -R "~{ref_fasta}" \
            -V "~{input_vcf}" \
            -sn "~{sample_names_file}" \
            -L "~{contig}" \
            -O "~{output_vcf_name}"
    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.6.1.0"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: memory_mb + " MiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File output_vcf = output_vcf_name
        File output_vcf_index = output_vcf_name + ".tbi"
    }
}

# Combine the per-contig VCFs back into a single file.
task ConcatenateVcfs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indexes

        String output_basename

        Int disk_size_gb = ceil(2 * size(input_vcfs, "GiB") + 10)
        Int cpu = 1
        Int memory_mb = 6000
        Int preemptible = 0
    }

    String output_vcf_name = output_basename + ".vcf.gz"

    command <<<
        set -euo pipefail

        bcftools concat --naive -Oz -o "~{output_vcf_name}" ~{sep=' ' input_vcfs}
        bcftools index -t ~{output_vcf_name}
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: memory_mb + " MiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File output_vcf = output_vcf_name
        File output_vcf_index = output_vcf_name + ".tbi"
    }
}
