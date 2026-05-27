# AGENTS.md

## What this repo is
- `terra-scientific-pipelines-service` ("Teaspoons") is a Spring Boot service that runs Terra scientific pipelines users cannot run directly (often due to protected reference data).
- Gradle multi-module layout (`settings.gradle`): `service` (runtime API), `client` (generated Java client), `rawls-client` (generated Rawls SDK), `python-client` (generated thin Python client), plus scripts.
- API contract lives in `common/openapi.yml`; server/client code is generated from it (`service/generators.gradle`, `client/java-client.gradle`, `python-client/teaspoons-client.gradle`).

## Runtime architecture and data flow
- Entry point: `service/src/main/java/bio/terra/pipelines/App.java`; component scan pulls in TCL (Stairway, IAM, tracing, migrate) and `bio.terra.pipelines`.
- Pipeline execution is intentionally two-phase: `prepare` then `start` (`PipelineRunsApiController` + `PipelineRunsService`).
- `prepare` validates inputs/quota, checks cloud input access, writes `pipeline_runs` + `pipeline_inputs`, and may return signed upload URLs.
- `start` requires `PREPARING` state, then submits a Stairway flight (`RunWdlBasedPipelineJobFlight` v20260428) that orchestrates Rawls/Cromwell steps (quota run -> input QC -> main workflow -> outputs).
- Success path writes outputs + quota to DB; failure path is handled via Stairway hooks/mark-failed logic (`service/IMPLEMENTATION_NOTES.md`).
- Output delivery is a separate Stairway flight (`DeliverDataToGcsFlight` v20260409): create delivery record -> copy outputs -> mark success -> best-effort source cleanup.

## Persistent model and boundaries
- DB schema is Liquibase-managed from `service/src/main/resources/db/changelog.xml`; add new changesets, do not rewrite historical ones.
- `PipelineRun` (`service/src/main/java/bio/terra/pipelines/db/entities/PipelineRun.java`) is source of truth for user-visible run status; Stairway job metadata can expire.

## External integrations you must preserve
- Sam (`service/src/main/java/bio/terra/pipelines/dependencies/sam/SamService.java`): authn/authz, admin checks, user pet SA tokens, proxy group lookup.
- Rawls (`service/src/main/java/bio/terra/pipelines/dependencies/rawls/RawlsService.java`): workspace metadata, method config validation/update, entity writes, submission status polling.
- GCS (`service/src/main/java/bio/terra/pipelines/dependencies/gcs/GcsService.java`): signed URL generation, IAM permission probes, copy/delete, requester-pays detection.
- Notifications (`service/src/main/java/bio/terra/pipelines/notifications/NotificationService.java`): publishes Pub/Sub messages for Thurloe-backed emails using templates in `notification-templates/`.

## Repo-specific coding conventions
- Flight versioning is strict: breaking Flight changes require new dated package `vYYYYMMDD`; keep old versions until no flights are in progress (`service/src/main/java/bio/terra/pipelines/stairway/flights/README.md`).
- New pipelines must be wired in multiple places: `PipelinesEnum`, `pipelines-config.yml`, DB seed/config rows, and `PipelineRunsService.startPipelineRun` switch.
- Keep generated sources out of manual edits (`service/build/swagger-code`, `client/build/swagger-code`, `python-client/generated`). Edit `common/openapi.yml` or generator configs instead.

## Making changes
- Any change that affects API contract, pipeline definition schema, or runtime behavior must be accompanied by unit/integration tests that validate the new behavior and guard against regressions. We require 80%+ coverage for new code paths.
