"""
Deploy updates for the Imputation marketing UI

Recursively uploads files from the current directory to the specified GCS bucket
based on the target environment. Certain files are excluded from the upload
(e.g. .DS_Store, deploy.py, readme.md).

Destination URLs:
* Development: https://allofus-anvil-imputation.dsde-dev.broadinstitute.org
* Production: https://allofus-anvil-imputation.terra.bio

See readme.md for more context.

Usage:
    python3 deploy.py --environment dev
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

# Define destination configurations
ENV_CONFIG = {
    "dev": {
        "bucket": "gs://allofus-anvil-imputation-dev",
        "project": "broad-dsde-dev"
    },
    "prod": {
        "bucket": "gs://allofus-anvil-imputation-prod",
        "project": "broad-dsde-prod"
    }
}

EXCLUDED_FILES = {".DS_Store", "deploy.py", "readme.md", "README.md"}

def validate_working_directory():
    expected_suffix = "ui/allofus-anvil-imputation"
    cwd = os.getcwd()
    if not cwd.endswith(expected_suffix):
        print(f"Error: This script must be run from a directory ending in '{expected_suffix}'")
        print(f"Current directory: {cwd}")
        sys.exit(1)

def get_current_gcloud_account() -> str:
    try:
        result = subprocess.run(
            ["gcloud", "auth", "list", "--filter=status:ACTIVE", "--format=value(account)"],
            check=True,
            stdout=subprocess.PIPE,
            text=True
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        print("Error: Failed to retrieve active gcloud account.")
        sys.exit(1)

def check_account_for_prod():
    account = get_current_gcloud_account()
    if not account.endswith("@firecloud.org"):
        print("Error: For 'prod' deployments, the active gcloud account must end in '@firecloud.org'.")
        print(f"Current account: {account}")
        sys.exit(1)

def should_exclude(path: str) -> bool:
    return os.path.basename(path) in EXCLUDED_FILES

def copy_files_to_temp_dir(source_dir=".") -> str:
    temp_dir = tempfile.mkdtemp()
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            full_path = os.path.join(root, file)
            rel_path = os.path.relpath(full_path, source_dir)
            if should_exclude(rel_path):
                continue
            dest_path = os.path.join(temp_dir, rel_path)
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)
            shutil.copy2(full_path, dest_path)
    return temp_dir

def run_gcloud_copy(environment: str):
    if environment not in ENV_CONFIG:
        print(f"Error: Unknown environment '{environment}'. Available: {', '.join(ENV_CONFIG.keys())}")
        sys.exit(1)

    if environment == "prod":
        check_account_for_prod()

    validate_working_directory()

    bucket = ENV_CONFIG[environment]["bucket"]
    project = ENV_CONFIG[environment]["project"]

    print("Collecting files to upload (excluding ignored files)...")
    temp_dir = copy_files_to_temp_dir()

    command = [
        "gcloud", "storage", "cp", "-r", ".", bucket,
        "--project", project,
         "--cache-control", "max-age=600, must-revalidate"
    ]

    try:
        print(f"Uploading from temporary directory: {temp_dir}")
        subprocess.run(command, check=True, cwd=temp_dir)
        print("✅ Upload completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"❌ Upload failed (code {e.returncode})")
        sys.exit(e.returncode)
    finally:
        shutil.rmtree(temp_dir)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-e", "--environment", type=str, required=True, help="Target environment (e.g. dev, prod)")
    args = parser.parse_args()
    run_gcloud_copy(args.environment)

if __name__ == "__main__":
    main()
