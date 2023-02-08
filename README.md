# Terra Scientific Pipelines Service

[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=DataBiosphere_terra-scientific-pipelines-service&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=DataBiosphere_terra-scientific-pipelines-service)

## Overview

Terra Scientific Pipelines Service, or teaspoons (tsps), facilitates running a number of defined scientific pipelines 
on behalf of users that users can't run themselves in Terra. The most common reason for this is that the pipeline 
accesses proprietary data that users are not allowed to access directly, but that may be used as e.g. a reference panel 
for imputation.

## Supported pipelines

Current supported pipelines are:
- [in development] Imputation (TODO add link/info)

## Architecture

TODO

## Development

This codebase is in initial development.

To run locally:
1. Make sure you have the requirements (below) installed. We recommend IntelliJ as an IDE.
2. Clone the repo
3. Run the commands in `scripts/postgres-init.sql` in your local postgres instance. You will need to be authenticated to access Vault.
4. Run `scripts/write-config.sh`
5. Run `./gradlew bootRun`
6. Navigate to [http://localhost:8080/#](http://localhost:8080/#)

### Requirements

This service is written in Java 17, and uses Postgres 13. 

To run locally, you'll also need: 
- jq - install with `brew install jq`
- vault - see DSP's setup instructions [here](https://docs.google.com/document/d/11pZE-GqeZFeSOG0UpGg_xyTDQpgBRfr0MLxpxvvQgEw/edit#heading=h.1k9ij99wmle2)
  - Note that for Step 7, "Create a GitHub Personal Access Token", you'll want to choose the "Tokens (classic)" option, not the fine-grained access token option.

### Tech stack

- Java 17 temurin
- Postgres 13.1
- Gradle - build automation tool
- SonarQube - static code security and coverage
- Trivy - security scanner for docker images
- Jib - docker image builder for Java