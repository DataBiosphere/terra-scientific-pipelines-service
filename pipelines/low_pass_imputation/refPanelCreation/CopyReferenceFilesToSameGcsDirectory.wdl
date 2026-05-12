version 1.0

workflow CopyReferenceFilesToSameDirectory {
    input {
        Array[String] contigs
        Array[File] sites_vcfs
        Array[File] sites_vcf_indices
        Array[File] sites_tables
        Array[File] sites_table_indices
        File reference_chunks
        String google_cloud_path
    }

    scatter (index in range(length(contigs))) {
        call SplitReferenceChunksAndCopy {
            input:
                contig = contigs[index],
                reference_chunks = reference_chunks,
                google_cloud_path = google_cloud_path
        }

        call RenameReferenceFilesAndCopy {
            input:
                contig = contigs[index],
                sites_vcf = sites_vcfs[index],
                sites_vcf_index = sites_vcf_indices[index],
                sites_table = sites_tables[index],
                sites_table_index = sites_table_indices[index],
                google_cloud_path = google_cloud_path
        }

    }

    output {
        Array[File] SplitReferenceChunksAndCopy_stdout = SplitReferenceChunksAndCopy.stdout
        Array[File] SplitReferenceChunksAndCopy_stderr = SplitReferenceChunksAndCopy.stderr
        Array[File] RenameReferenceFilesAndCopy_stdout = RenameReferenceFilesAndCopy.stdout
        Array[File] RenameReferenceFilesAndCopy_stderr = RenameReferenceFilesAndCopy.stderr
    }
}

task SplitReferenceChunksAndCopy {
    input {
        String contig
        File reference_chunks
        String google_cloud_path
        Int disk_size = ceil(size(reference_chunks, "GiB") + 10)
    }

    command <<<
        set -euo pipefail

        grep "_~{contig}_" ~{reference_chunks} > reference_chunks.~{contig}.txt
        gcloud storage cp reference_chunks.~{contig}.txt ~{google_cloud_path}
    >>>

    runtime {
        docker: "google/cloud-sdk:562.0.0-slim"
        memory: "4 GiB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

    output{
        File stdout = stdout()
        File stderr = stderr()
    }
}

task RenameReferenceFilesAndCopy {
    input {
        String contig
        File sites_vcf
        File sites_vcf_index
        File sites_table
        File sites_table_index

        String google_cloud_path
        Int disk_size = ceil(2*(size(sites_vcf, "GiB") + size(sites_vcf_index, "GiB") + size(sites_table, "GiB") + size(sites_table_index, "GiB"))  + 10)
    }

    command <<<
        set -euo pipefail

        # rename all the files
        mv ~{sites_vcf} sites.~{contig}.vcf.gz
        mv ~{sites_vcf_index} sites.~{contig}.vcf.gz.tbi
        mv ~{sites_table} sites_table.~{contig}.gz
        mv ~{sites_table_index} sites_table.~{contig}.gz.tbi

        # copy files to google cloud path
        gcloud storage cp site* ~{google_cloud_path}
    >>>

    runtime {
        docker: "google/cloud-sdk:562.0.0-slim"
        memory: "4 GiB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }

    output{
        File stdout = stdout()
        File stderr = stderr()
    }
}
