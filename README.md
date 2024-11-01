# Terra Scientific Pipelines Service

[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=DataBiosphere_terra-scientific-pipelines-service&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=DataBiosphere_terra-scientific-pipelines-service)

## Overview

Terra Scientific Pipelines Service, or teaspoons, facilitates running a number of defined scientific pipelines 
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
7. If this is your first time deploying to any environment, be sure to use the admin endpoint `/api/admin/v1/updatePipelineWorkspaceId/{pipelineName}/{workspaceId}` to set your pipeline's workspace id. Workspace id can be found through the terra ui workspace dashboard or through the Rawls GET workspace endpoint.


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

### (Optional) Install pre-commit hooks
1. [scripts/git-hooks/pre-commit] has been provided to help ensure all submitted changes are formatted correctly.  To install all hooks in [scripts/git-hooks], run:
```bash
git config core.hooksPath scripts/git-hooks
```

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

To connect to the Teaspoons database, we have a script in [dsp-scripts](https://github.com/broadinstitute/dsp-scripts) that 
does all the setup for you. Clone that repo and make sure you're either on Broad Internal wifi or connected
to the VPN. Then run the following command:

```shell
./db/psql-connect.sh dev tsps
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

### Running the end-to-end tests

The end-to-end test is specified in `.github/workflows/run-e2e-tests.yaml`. It calls [the test script defined 
in the dsp-reusable-workflows repo](https://github.com/broadinstitute/dsp-reusable-workflows/blob/main/e2e-test/tsps_e2e_test.py).

The end-to-end test is automatically run nightly on the dev environment. 

To run the test against a specific feature branch:
1. Grab the image tag for your feature branch. 
> If you've opened a PR, you can find the image tag as follows:
> - go to the Bump, Tag, Publish, and Deploy workflow that's triggered each time you push to your branch
> - From there, go to the tag-publish-docker-deploy task
> - Expand the "Construct docker image name and tag" step
> - The first line should contain the image tag, something like "0.0.81-6761487".
2. Navigate to the [e2e-test GHA workflow](https://github.com/DataBiosphere/terra-scientific-pipelines-service/actions/workflows/run-e2e-tests.yaml)
3. Click on the "Run workflow" button and select your branch from the dropdown
 - Enter the image tag from step 1 in the "Custom image tag" field
 - If you've updated the end-to-end test in the dsp-resuable-workflows repo, enter either a commit hash or your git 
branch name. If you don't need to change the test, leave the default as main.
4. Click the green "Run workflow" button.

## Python clients
We publish a "thin", auto-generated Python client that wraps the Teaspoons APIs. This client is published to 
[PyPi](https://pypi.org/project/terra-scientific-pipelines-service-api-client/) and can be installed with 
`pip install teaspoons_client`, although this is not meant to be user-facing. The thin api client is generated from 
the OpenAPI spec in the `openapi` directory. 

Publishing occurs automatically when a new version of the service is deployed, via the 
[release-python-client GHA](https://github.com/DataBiosphere/terra-scientific-pipelines-service/blob/main/.github/workflows/release-python-client.yml). 

We also have a user-facing, "thick" CLI whose code lives in a separate repository: [DataBiosphere/terra-scientific-pipelines-service-cli](https://github.com/DataBiosphere/terra-scientific-pipelines-service-cli). 
