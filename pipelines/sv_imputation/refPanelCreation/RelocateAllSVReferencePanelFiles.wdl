version 1.0

import "RelocateJsonFiles.wdl" as RelocateJsonFiles

# Relocates all of the final SV reference panel outputs to their destinations: every
# gs:// file referenced by the chunked_panel_json and pop_panel_resources_json manifests
# (via the shared RelocateFiles task), both manifests themselves, and the
# preprocess_panel_bubble_split_sites_only_vcf file + index produced separately by
# CreatePanelAuxiliaryFiles.
workflow RelocateAllSVReferencePanelFiles {
    input {
        File chunked_panel_json
        File pop_panel_resources_json

        String preprocess_panel_bubble_split_sites_only_vcf
        String preprocess_panel_bubble_split_sites_only_vcf_idx

        String auxiliary_files_destination_gcs_path
        String json_destination_gcs_path
        String panel_bin_files_destination_gcs_path

        Boolean move_files = false
        Boolean dry_run = false
    }

    call RelocateJsonFiles.RelocateFiles as RelocateChunkedPanelJsonManifest {
        input:
            input_json = chunked_panel_json,
            destination_gcs_path = panel_bin_files_destination_gcs_path,
            json_destination_gcs_path = json_destination_gcs_path,
            move_files = move_files,
            dry_run = dry_run
    }

    call RelocateJsonFiles.RelocateFiles as RelocatePopPanelResourcesJsonManifest {
        input:
            input_json = pop_panel_resources_json,
            destination_gcs_path = auxiliary_files_destination_gcs_path,
            json_destination_gcs_path = json_destination_gcs_path,
            move_files = move_files,
            dry_run = dry_run
    }

    call RelocateFilePair as RelocateSitesOnlyVcf {
        input:
            source_file = preprocess_panel_bubble_split_sites_only_vcf,
            source_file_index = preprocess_panel_bubble_split_sites_only_vcf_idx,
            destination_gcs_path = auxiliary_files_destination_gcs_path,
            move_files = move_files,
            dry_run = dry_run
    }

    output {
        File chunked_panel_updated_json = RelocateChunkedPanelJsonManifest.updated_json
        File chunked_panel_original_paths_fofn = RelocateChunkedPanelJsonManifest.original_paths_fofn
        String chunked_panel_updated_json_gcs_path = RelocateChunkedPanelJsonManifest.updated_json_gcs_path
        File pop_panel_resources_updated_json = RelocatePopPanelResourcesJsonManifest.updated_json
        File pop_panel_resources_original_paths_fofn = RelocatePopPanelResourcesJsonManifest.original_paths_fofn
        String pop_panel_resources_updated_json_gcs_path = RelocatePopPanelResourcesJsonManifest.updated_json_gcs_path
        String relocated_preprocess_panel_bubble_split_sites_only_vcf = RelocateSitesOnlyVcf.relocated_file
        String relocated_preprocess_panel_bubble_split_sites_only_vcf_idx = RelocateSitesOnlyVcf.relocated_file_index
    }
}

# Copies or moves a file and its index to a gs:// destination (bucket and prefix).
# source_file/source_file_index are gs:// paths passed through as Strings (not File)
# so that gcloud storage mv/cp operates directly on the source and destination URIs -
# a true, efficient server-side move (deleting the real source object) rather than a
# download-then-upload of a localized copy.
task RelocateFilePair {
    input {
        String source_file
        String source_file_index
        String destination_gcs_path
        Boolean move_files
        Boolean dry_run

        Int disk_size_gb = 10
        Int cpu = 1
        Int memory_mb = 4000
        String cloud_sdk_docker = "google/cloud-sdk:562.0.0-slim"
    }

    # strip any trailing slash(es) so the destination paths are constructed cleanly
    String dest = sub(destination_gcs_path, "/+$", "")
    String action = if move_files then "mv" else "cp"
    String file_name = basename(source_file)
    String file_index_name = basename(source_file_index)

    command <<<
        set -euo pipefail

        if ~{if dry_run then "true" else "false"}; then
            echo "[dry run] would ~{action}: ~{source_file} -> ~{dest}/~{file_name}"
            echo "[dry run] would ~{action}: ~{source_file_index} -> ~{dest}/~{file_index_name}"
        else
            gcloud storage ~{action} ~{source_file} ~{dest}/~{file_name}
            gcloud storage ~{action} ~{source_file_index} ~{dest}/~{file_index_name}
        fi
    >>>

    runtime {
        docker: cloud_sdk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 0
    }

    output {
        String relocated_file = "~{dest}/~{file_name}"
        String relocated_file_index = "~{dest}/~{file_index_name}"
    }
}
