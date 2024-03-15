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

[WIP architecture doc](https://docs.google.com/document/d/1dAPwOG2z1h0B5CszeQ0DfyToniNV_3y1OBV7x7L8ofI/edit?usp=sharing)
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
- Java 17 - can be installed manually or through IntelliJ which will do it for you when importing the project
- Postgres 13 - multiple solutions here as long as you have a postgres instance running on localhost:5432 the local app will connect appropriately
  - Postgres.app https://postgresapp.com/
  - Brew https://formulae.brew.sh/formula/postgresql@13

### Tech stack

- Java 17 temurin
- Postgres 13.1
- Gradle - build automation tool
- SonarQube - static code security and coverage
- Trivy - security scanner for docker images
- Jib - docker image builder for Java

### Local development

To run locally:
1. Make sure you have the requirements installed from above. We recommend IntelliJ as an IDE.
2. Clone the repo (if you see broken inputs build the project to get the generated sources)
3. Run the commands in `scripts/postgres-init.sql` in your local postgres instance. You will need to be authenticated to access Vault.
4. Run `scripts/write-config.sh`
5. Run `./gradlew bootRun` to spin up the server.
6. Navigate to [http://localhost:8080/#](http://localhost:8080/#)

#### Local development with debugging
If using Intellij (only IDE we use on the team), you can run the server with a debugger. Follow
the steps above but instead of running `./gradlew bootRun` to spin up the server, you can run
(debug) the App.java class through intellij and set breakpoints in the code.  Be sure to set the
GOOGLE_APPLICATION_CREDENTIALS=config/tsps-sa.json in the Run/Debug configuration Environment Variables.

### Running Tests/Linter Locally
- Testing
  - Run `./gradlew service:test` to run tests
- Linting
  - Run `./gradlew spotlessCheck` to run linter checks 
  - Run `./gradlew :service:spotlessApply` to apply fix any issues the linter finds

### Running SonarQube locally

[SonarQube](https://www.sonarqube.org) is a static analysis code that scans code for a wide
range of issues, including maintainability and possible bugs. Get more information from
[DSP SonarQube Docs](https://dsp-security.broadinstitute.org/appsec-team-internal/appsec-team-internal/security-activities/sast-1#)

If you get a build failure due to
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

### Connecting to the database

To connect to the TSPS database, we have a script in [dsp-scripts](https://github.com/broadinstitute/dsp-scripts) that 
does all the setup for you. Clone that repo and make sure you're either on Broad Internal wifi or connected
to the VPN. Then run the following command:

```shell
./firecloud/psql-connect.sh dev tsps
```

### Deploying to dev

Upon merging to main, the dev environment will be automatically deployed via the GitHub Action [Bump, Tag, Publish, and Deploy](https://github.com/DataBiosphere/terra-scientific-pipelines-service/actions/workflows/tag-publish.yml)
(that workflow is defined [here](https://github.com/DataBiosphere/terra-scientific-pipelines-service/blob/main/.github/workflows/tag-publish.yml)). 

The two tasks `report-to-sherlock` and `set-version-in-dev` will prompt Sherlock to deploy the new version to dev. 
You can check the status of the deployment in [Beehive](https://beehive.dsp-devops.broadinstitute.org/apps/tsps) and in 
[ArgoCD](https://ap-argocd.dsp-devops.broadinstitute.org/applications/ap-argocd/tsps-dev).

For more information about deployment to dev, check out DevOps' [excellent documentation](https://docs.google.com/document/d/1lkUkN2KOpHKWufaqw_RIE7EN3vN4G2xMnYBU83gi8VA/).

### Tracing

We use [OpenTelemetry](https://opentelemetry.io/) for tracing, so that every request has a tracing span that can 
be viewed in [Google Cloud Trace](https://cloud.google.com/trace). (This is not yet fully set up here - to be done in TSPS-107). 
See [this DSP blog post](https://broadworkbench.atlassian.net/wiki/x/AoGlrg) for more info.
