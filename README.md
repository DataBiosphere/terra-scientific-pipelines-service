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

### Requirements

This service is written in Java 17, and uses Postgres 13. 

### Tech stack

- Java 17 temurin
- Postgres 13.1
- Gradle - build automation tool
- SonarQube - static code security and coverage
- Trivy - security scanner for docker images
- Jib - docker image builder for Java