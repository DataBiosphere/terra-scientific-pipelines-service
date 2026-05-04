#!/usr/bin/env bash
# usage from root repo directory: ./scripts/local-dev/copy_web_config.sh
# this script copies the LocalWebConfig.java from this directory to next to App.java, where it will be ignored by git
# and will allow local development (avoid CORS issues when running the service and UI locally)

cp scripts/local-dev/LocalWebConfig.java service/src/main/java/bio/terra/pipelines/LocalWebConfig.java
