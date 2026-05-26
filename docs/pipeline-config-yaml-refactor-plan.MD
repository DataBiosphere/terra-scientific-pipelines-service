# Pipeline Config YAML Refactor Plan

This document is the working plan for migrating pipeline configuration from the database to `pipelines-config.yml`, while introducing an internal canonical `pipeline_key`.

## Summary
We will move most pipeline definition metadata out of the database and into `service/src/main/resources/pipelines-config.yml`, and migrate internal runtime/persistence lookups from numeric `pipeline_id` to string `pipeline_key`.

Canonical key format:

`pipeline_key = {pipelineName}_v{pipelineVersion}`
- all lowercase
- underscore as separator between pipeline name (which itself has underscores) and version

Example: `array_imputation_v3`

## Scope
### In scope
- Move full contents of:
  - `pipeline_input_definitions`
  - `pipeline_output_definitions`
  - `pipeline_quotas`
  - all pipeline-definition fields currently stored in `pipelines` except the following, which should be copied into a new `pipeline_runtime_metadata` table for persistence of runtime-specific data:
    - `pipeline_id` (legacy numeric ID, for backward compatibility)
    - `pipeline_key` (new string ID composed of existing name and version fields: `{name}_v{version}`)
    - `tool_version`
    - `workspace_billing_project`
    - `workspace_name`
    - `workspace_google_project`
    - `workspace_storage_container_name`
    - `hidden`
  - the admin updatePipeline endpoint already only touches the fields above, so no API contract changes are needed.
- Introduce and propagate internal `pipeline_key` usage in service/runtime layers.
  - The `pipeline_runs` table will retain both `pipeline_id` and `pipeline_key` for the foreseeable future to support backward compatibility during migration.
- Add additive database migrations to support compatibility during rollout.
- Maintain backward compatibility for in-flight Stairway jobs using versioned flights where needed.

### Out of scope (explicit)
- No API endpoint changes.
- No API request/response field changes.
- No `common/openapi.yml` contract changes.
- No generated client regeneration or client migration work.

External API behavior remains based on `pipelineName` and `pipelineVersion`.

## Goals
- Make `pipelines-config.yml` the source of truth for pipeline definitions.
- Reduce DB-coupled config complexity and simplify pipeline definition updates.
- Keep runtime stable during migration.

## No Dual-Source Decision
- This refactor will use a **single release cutover** from database-backed pipeline definitions to YAML-backed pipeline definitions.
- There will be **no intermediate dual-source runtime mode** where some pods or code paths read definitions from the database while others read from YAML.
- External API behavior remains unchanged throughout the cutover; callers continue to use `pipelineName` and `pipelineVersion`.
- To avoid mixed-source behavior during deployment, release sequencing must ensure that old and new runtime definition sources are not serving traffic at the same time.
- Rollback is a **release rollback decision**, not a runtime source toggle. If cutover fails, the rollback procedure must restore the previous application version and any required database state for legacy definition reads.

## Non-goals
- Providing rollback controls via feature flags.
- Redesigning quota semantics.
- Reworking user-visible identifiers in API payloads.
- Unrelated refactors across Sam/Rawls/GCS integration boundaries.

## Target end state
- Pipeline definition metadata (inputs, outputs, quotas/defaults, and static pipeline fields) is loaded from YAML at startup.
- Internal service logic resolves pipeline definitions by `pipeline_key`.
- Legacy DB config tables/columns are no longer required for active read paths.

## Phased plan

### Phase 0: Design guardrails and invariants
- Record constraints in a design note/architecture decision record (ADR):
  - API contract unchanged.
  - Additive DB migration only.
  - Backward-compatible Stairway strategy.
- Define `pipeline_key` normalization rules (case, separator, version formatting).

### Phase 1: YAML schema and typed loader
- Extend `service/src/main/resources/pipelines-config.yml` with canonical sections for:
  - pipeline metadata
  - input definitions
  - output definitions
  - quota/default config currently tied to definition tables
- Add/extend typed config classes and validation at startup.
- Unit tests validate config (check for duplicate keys, missing required fields, bad references).

### Phase 2: Additive DB compatibility migration
- Add new Liquibase changesets in `service/src/main/resources/db/changelog.xml` to:
  - add `pipeline_key` columns where runtime/persistence needs identity linkage
  - backfill from current `(pipeline_name, pipeline_version)` values
  - add indexes/uniqueness constraints for efficient key lookup
- Keep legacy columns/tables in place.

### Phase 3: Pre-cutover internal refactor and readiness
- Update internal services (prepare/start/quota/input-output) to work with key-centric model and YAML-backed pipeline definitions.
- Keep controller/API mapping unchanged.
- Persist all execution-critical runtime metadata required after prepare/start in `pipeline_runtime_metadata` and `pipeline_runs`, so execution does not depend on removed DB definition rows.
- Add unit tests that fail if:
  - YAML contains duplicate or malformed `pipeline_key` values
  - required pipeline metadata, input definitions, output definitions, or quota entries are missing
- Add runtime validation at application startup to fail fast if:
  - runtime metadata cannot be correlated to configured pipelines
- Ensure that:
  - all pipelines present in production DB definitions are represented in YAML
  - `pipeline_key` backfill is complete and unique
  - historic and active `pipeline_runs` rows needed by the service remain resolvable after cutover

### Phase 4: Single-release cutover and Stairway compatibility
- For breaking flight input-map changes, introduce a new dated flight package per repo convention.
- For breaking step working map changes, version step classes.
- Keep old flight versions active until no in-flight jobs depend on them.
- Execute cutover as a coordinated single release:
  - complete DB migrations and backfills before enabling the new application version
  - verify no incompatible in-flight Stairway jobs remain
  - deploy the application version that reads pipeline definitions only from YAML
  - run immediate post-deploy health checks before considering cutover complete

### Phase 5: Post-cutover validation and rollback stance
- Validate cutover success with targeted checks:
  - application startup succeeds with minimal YAML validation enabled
  - pipeline list/details APIs return expected content without relying on DB definition tables
  - prepare/start/result flows succeed for representative pipelines
  - admin update flows continue to update `pipeline_runtime_metadata` correctly
- Treat rollback as a release-management operation, not a source toggle:
  - if cutover fails, revert to the prior application version
  - restore any required database state for legacy definition reads if the failed release made incompatible cleanup changes
  - defer destructive cleanup until the YAML-backed release has completed a soak period

### Phase 6: Cleanup
- Remove dead read paths from definition tables.
- Add final Liquibase cleanup changesets to drop obsolete config tables/columns after soak period.

## Testing strategy
- Unit tests:
  - `pipeline_key` generation and parsing
  - YAML binding/validation
  - adding new pipelines for testing
- Migration tests:
  - Liquibase backfill correctness
  - uniqueness/index constraints
- Stairway tests:
  - old/new flight compatibility
  - working-map field compatibility
- Regression/smoke tests:
  - unchanged API payloads and endpoint behavior

## Risks and mitigations
- In-flight jobs expect old identifiers
  - Mitigation: versioned flights during transition.
- Mixed-source behavior during deployment causes inconsistent pipeline resolution
  - Mitigation: coordinated single-release cutover with no simultaneous serving of DB-backed and YAML-backed definition readers.
- Cutover fails and rollback is unclear without dual-source support
  - Mitigation: define and rehearse a release rollback runbook before production deployment.
- Backfill or YAML omissions leave runs or pipeline metadata unresolvable after cutover
  - Mitigation: require explicit pre-cutover parity checks and startup fail-fast validation.

## Incremental PR sequence
1. Design doc + `pipeline_key` rules + YAML schema scaffolding.
2. Additive DB migrations for `pipeline_key` + backfill + indexes.
3. Internal refactor for key-centric resolution, persisted runtime metadata, and startup validation.
4. Stairway compatibility updates (flight and step updates).
5. YAML source-of-truth implementation and single-release cutover preparation.
6. Production cutover release and post-cutover validation.
7. Legacy DB definition path removal + schema cleanup.

---

Owner: Teaspoons service team
Status: Draft for review and iteration

