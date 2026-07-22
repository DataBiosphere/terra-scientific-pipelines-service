version 1.0

# Shared task to copy a file and its index to a gs:// destination (bucket and prefix).
# Imported by the wdls in this directory to optionally delocalize their outputs.
task CopyToCloud {
    input {
        File source_file
        File source_file_index
        String copy_to_cloud_dest

        Int disk_size_gb = ceil((size(source_file, "GiB") + size(source_file_index, "GiB"))) + 10
        Int cpu = 1
        Int memory_mb = 4000
        String cloud_sdk_docker = "google/cloud-sdk:562.0.0-slim"
    }

    # strip any trailing slash(es) so the destination paths are constructed cleanly
    String dest = sub(copy_to_cloud_dest, "/+$", "")
    String file_name = basename(source_file)
    String file_index_name = basename(source_file_index)

    command <<<
        set -euo pipefail

        gcloud storage cp ~{source_file} ~{dest}/~{file_name}
        gcloud storage cp ~{source_file_index} ~{dest}/~{file_index_name}
    >>>

    runtime {
        docker: cloud_sdk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 0
    }

    output {
        File copied_file = "~{dest}/~{file_name}"
        File copied_file_index = "~{dest}/~{file_index_name}"
    }
}
