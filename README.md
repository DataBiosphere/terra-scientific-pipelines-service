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

[WIP architecture doc](https://docs.google.com/document/d/1UfkWpBInSqnd5gDF-8KuR6pPC5Nsih4slqeI1HK63Zg/edit?usp=sharing)
[Linked LucidChart](https://lucid.app/lucidchart/2f067b5e-2d40-41b4-a5f3-a9dc72d83820/edit?viewport_loc=-72%2C25%2C1933%2C1133%2C0_0&invitationId=inv_97522cca-1b6d-44fe-9552-8f959d410dd7)

## Development

This codebase is in initial development.

### Requirements

This service is written in Java 17, and uses Postgres 13.

To run locally, you'll also need:
- jq - install with `brew install jq`
- vault - see DSP's setup instructions [here](https://docs.google.com/document/d/11pZE-GqeZFeSOG0UpGg_xyTDQpgBRfr0MLxpxvvQgEw/edit#heading=h.1k9ij99wmle2)
  - Note that for Step 7, "Create a GitHub Personal Access Token", you'll want to choose
    the "Tokens (classic)" option, not the fine-grained access token option.

### Tech stack

- Java 17 temurin
- Postgres 13.1
- Gradle - build automation tool
- SonarQube - static code security and coverage
- Trivy - security scanner for docker images
- Jib - docker image builder for Java

### Local development

To run locally:
1. Make sure you have the requirements (below) installed. We recommend IntelliJ as an IDE.
2. Clone the repo
3. Run the commands in `scripts/postgres-init.sql` in your local postgres instance. You will need to be authenticated to access Vault.
4. Run `scripts/write-config.sh`
5. Run `./gradlew bootRun`
6. Navigate to [http://localhost:8080/#](http://localhost:8080/#)

### Running SonarQube locally

[SonarQube](https://www.sonarqube.org) is a static analysis code that scans code for a wide
range of issues, including maintainability and possible bugs. If you get a build failure due to
SonarQube and want to debug the problem locally, you need to get the sonar token from vault
before running the gradle task.

```shell
export SONAR_TOKEN=$(vault read -field=sonar_token secret/secops/ci/sonarcloud/tsps)
./gradlew sonarqube
```

Running this task produces no output unless your project has errors. To
generate a report, run using `--info`:

```shell
./gradlew sonarqube --info
```