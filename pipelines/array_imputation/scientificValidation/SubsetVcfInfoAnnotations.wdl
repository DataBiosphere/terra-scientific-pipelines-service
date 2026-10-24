version 1.0

workflow SubsetVcfInfoAnnotations {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indexes
        String output_basename
        Int preemptible = 0
    }

    if (length(input_vcfs) != length(input_vcf_indexes)) {
        call FailOnMismatchedInputs {
            input:
                output_basename = output_basename,
                preemptible = preemptible
        }
    }

    scatter (i in range(length(input_vcfs))) {
        call KeepAllowedInfoAnnotations {
            input:
                input_vcf = input_vcfs[i],
                input_vcf_index = input_vcf_indexes[i],
                output_basename = output_basename + "." + basename(input_vcfs[i], ".vcf.gz"),
                preemptible = preemptible
        }
    }

    call ConcatenateVcfs {
        input:
            input_vcfs = KeepAllowedInfoAnnotations.output_vcf,
            input_vcf_indexes = KeepAllowedInfoAnnotations.output_vcf_index,
            output_basename = output_basename,
            preemptible = preemptible
    }

    output {
        File merged_vcf = ConcatenateVcfs.output_vcf
        File merged_vcf_index = ConcatenateVcfs.output_vcf_index
    }
}

# Keep only the requested INFO annotations and drop everything else.
# This preserves the VCF structure while removing all non-whitelisted INFO fields.
task KeepAllowedInfoAnnotations {
    input {
        File input_vcf
        File input_vcf_index
        String output_basename

        Int disk_size_gb = ceil(2 * size(input_vcf, "GiB") + 20)
        Int cpu = 1
        Int memory_mb = 6000
        Int preemptible = 0
    }

    String output_vcf_name = output_basename + ".filtered.vcf.gz"

    command <<<
        set -euo pipefail

        all_info_ids=$(bcftools view -h "~{input_vcf}" | awk -F'[=,]' '/^##INFO=<ID=/{print $3}' | sort -u | paste -sd, -)
        remove_info_ids=$(printf '%s\n' "$all_info_ids" | tr ',' '\n' | grep -vE -x 'AF_afr|AF_ami|AF_amr|AF_asj|AF_eas|AF_fin|AF_mid|AF_nfe|AF_raw|AF_remaining|AF_sas' | paste -sd, -)

        if [ -n "$remove_info_ids" ]; then
            remove_info_id_list=$(printf '%s\n' "$remove_info_ids" | tr ',' '\n' | sed 's/^/INFO\//g' | paste -sd, -)
            bcftools annotate \
                -x "$remove_info_id_list" \
                -Oz \
                --write-index=tbi \
                -o "~{output_vcf_name}" \
                "~{input_vcf}"
        else
            bcftools view -Oz --write-index=tbi -o "~{output_vcf_name}" "~{input_vcf}"
        fi
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

# Combine the filtered per-file VCFs back into a single merged file.
task ConcatenateVcfs {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indexes

        String output_basename

        Int disk_size_gb = ceil(2.4 * size(input_vcfs, "GiB") + 10)
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

# Guard against mismatched list lengths. This is a runtime failure, not a WDL validation issue.
task FailOnMismatchedInputs {
    input {
        String output_basename
        Int preemptible = 0
    }

    command <<<
        echo "The VCF and VCF index arrays must be the same length." >&2
        exit 1
    >>>

    runtime {
        docker: "ubuntu:22.04"
        preemptible: preemptible
    }
}
