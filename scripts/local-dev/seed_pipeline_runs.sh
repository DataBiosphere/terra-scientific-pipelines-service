#!/usr/bin/env bash
set -euo pipefail

# Usage: ./trigger_pipeline_runs.sh <N>
if [ $# -ne 1 ]; then
  echo "Usage: $0 <number_of_runs>"
  exit 1
fi

N=$1
ACCESS_TOKEN=$(gcloud auth print-access-token)
SQL_FILE="random_status_updates.sql"

# Define valid statuses
STATUSES=("PREPARING" "SUCCEEDED" "FAILED" "RUNNING")

# Initialize the SQL file
{
  echo "-- Randomly assign statuses to pipeline runs"
  echo "-- Generated on $(date)"
  echo
} > "$SQL_FILE"

for i in $(seq 1 "$N"); do
  JOB_ID=$(uuidgen)

  # Pick a random status using Bash's $RANDOM
  RAND_INDEX=$((RANDOM % ${#STATUSES[@]}))
  STATUS=${STATUSES[$RAND_INDEX]}

  echo "Starting run $i (job ID: $JOB_ID, status: $STATUS)"

  # Call the API endpoint
  curl -s -o /dev/null -w "%{http_code}\n" -X POST \
    "http://localhost:8080/api/pipelineruns/v1/prepare" \
    -H "accept: */*" \
    -H "Authorization: Bearer $ACCESS_TOKEN" \
    -H "Content-Type: application/json" \
    -d "{
      \"jobId\": \"$JOB_ID\",
      \"pipelineName\": \"array_imputation\",
      \"pipelineVersion\": 1,
      \"pipelineInputs\": {
        \"multiSampleVcf\": \"test.vcf.gz\",
        \"outputBasename\": \"output\"
      },
      \"description\": \"Test run $i\",
      \"useResumableUploads\": false
    }" > /dev/null

  # Append SQL line
  echo "update pipeline_runs set status='${STATUS}' where job_id='${JOB_ID}';" >> "$SQL_FILE"

  echo " → Completed run $i"
done

echo
echo "All $N runs triggered."
echo "SQL file created: $SQL_FILE"
