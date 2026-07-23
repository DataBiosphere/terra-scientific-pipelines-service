version 1.0

import "CopyToCloud.wdl" as CopyToCloudTask

# This script is under review. It is not actively tested or maintained at this time.
workflow SplitMultiallelicsBcf {
    input {
        File input_bcf
        File input_bcf_index
        File ref_fasta
        File ref_fasta_index
        String contig
        String output_basename
        String? copy_to_cloud_dest
    }

    call SeparateMultiallelics {
        input:
            bcf_input = input_bcf,
            bcf_input_index = input_bcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            output_bcf_name = output_basename + ".multiallelic_split." + contig + ".bcf"
    }

    if (defined(copy_to_cloud_dest)) {
        call CopyToCloudTask.CopyToCloud as CopyToCloud {
            input:
                source_file = SeparateMultiallelics.output_bcf,
                source_file_index = SeparateMultiallelics.output_bcf_index,
                copy_to_cloud_dest = select_first([copy_to_cloud_dest])
        }
    }

    output {
        String multi_allelics_split_bcf = select_first([CopyToCloud.copied_file, SeparateMultiallelics.output_bcf])
        String multi_allelics_split_bcf_index = select_first([CopyToCloud.copied_file_index, SeparateMultiallelics.output_bcf_index])
    }
}

task SeparateMultiallelics {
    input {
        File bcf_input
        File bcf_input_index
        File ref_fasta
        File ref_fasta_index
        String output_bcf_name

        Int disk_size_gb =  ceil(3 * (size(bcf_input, "GiB") + size(bcf_input_index, "GiB"))) + 20
        Int cpu = 1
        Int memory_mb = 12000
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
    }


    command <<<
        set -e -o pipefail

        # split multiallelics w/left aligning and then recalculate AC, AN
        bcftools norm -m -f ~{ref_fasta} --rm-dup exact - ~{bcf_input} -Ou | bcftools +fill-tags - -Ob -o ~{output_bcf_name} -- -t AC,AN

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
