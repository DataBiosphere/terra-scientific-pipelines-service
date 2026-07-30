version 1.0

# Given an arbitrary JSON manifest containing gs:// file paths (nested at any depth,
# inside objects and/or arrays), copies or moves each referenced file to a single
# destination gs:// directory, and rewrites the manifest to point at the new locations.
#
# Files are relocated by basename only (destination_gcs_path/<basename>); if two
# distinct source paths would collide on the same destination basename, the task fails
# rather than silently overwriting one of them.
workflow RelocateJsonFiles {
    input {
        File input_json
        String destination_gcs_path
        Boolean move_files = false
        Boolean dry_run = false
    }

    call RelocateFiles {
        input:
            input_json = input_json,
            destination_gcs_path = destination_gcs_path,
            move_files = move_files,
            dry_run = dry_run
    }

    output {
        File updated_json = RelocateFiles.updated_json
        File original_paths_fofn = RelocateFiles.original_paths_fofn
    }
}

task RelocateFiles {
    input {
        File input_json
        String destination_gcs_path
        Boolean move_files
        Boolean dry_run

        Int cpu = 1
        Int memory_mb = 4000
        Int disk_size_gb = 10
        String cloud_sdk_docker = "google/cloud-sdk:562.0.0-slim"
    }

    # strip any trailing slash(es) so destination paths are constructed cleanly
    String dest = sub(destination_gcs_path, "/+$", "")
    String action = if move_files then "mv" else "cp"
    String dry_run_flag = if dry_run then "true" else "false"

    command <<<
        set -euo pipefail

        python3 <<CODE
        import json
        import re
        import subprocess

        with open("~{input_json}") as f:
            data = json.load(f)

        dest = "~{dest}"
        action = "~{action}"
        dry_run = "~{dry_run_flag}" == "true"
        gs_file_re = re.compile(r"^gs://.+[^/]$")

        # first pass: find every gs:// file path in the manifest (any nesting depth,
        # inside dicts and/or lists), and fail fast if two distinct source paths would
        # collide on the same destination basename
        basename_to_source = {}

        def collect(node):
            if isinstance(node, dict):
                for value in node.values():
                    collect(value)
            elif isinstance(node, list):
                for value in node:
                    collect(value)
            elif isinstance(node, str) and gs_file_re.match(node):
                basename = node.rsplit("/", 1)[-1]
                existing = basename_to_source.get(basename)
                if existing is not None and existing != node:
                    raise ValueError(
                        f"basename collision at destination: '{basename}' is shared by "
                        f"'{existing}' and '{node}'"
                    )
                basename_to_source[basename] = node

        collect(data)
        print(f"found {len(basename_to_source)} file(s) to {action}")

        # write out a FOFN of every original source path, before any relocation happens
        with open("original_paths.fofn", "w") as f:
            for source in sorted(basename_to_source.values()):
                f.write(source + "\n")

        # second pass: relocate each source file, building the old -> new path map
        path_map = {}
        for basename, source in sorted(basename_to_source.items()):
            new_path = f"{dest}/{basename}"
            path_map[source] = new_path
            if dry_run:
                print(f"[dry run] would {action}: {source} -> {new_path}")
                continue
            print(f"{action}: {source} -> {new_path}")
            subprocess.run(["gcloud", "storage", action, source, new_path], check=True)

        # third pass: rewrite the manifest with the new destination paths
        def relocate(node):
            if isinstance(node, dict):
                return {key: relocate(value) for key, value in node.items()}
            if isinstance(node, list):
                return [relocate(value) for value in node]
            if isinstance(node, str) and node in path_map:
                return path_map[node]
            return node

        with open("updated.json", "w") as f:
            json.dump(relocate(data), f, indent=2)
        CODE
    >>>

    runtime {
        docker: cloud_sdk_docker
        disks: "local-disk ${disk_size_gb} HDD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: 0
    }

    output {
        File updated_json = "updated.json"
        File original_paths_fofn = "original_paths.fofn"
    }
}
