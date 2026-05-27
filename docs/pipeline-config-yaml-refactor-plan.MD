# Pipeline Config YAML Refactor Plan

This plan migrates pipeline definitions from DB tables to `service/src/main/resources/pipelines-config.yml` and standardizes internal identity on `pipeline_key`.

## Summary
- Source of truth for pipeline definitions becomes YAML.
- Internal identity changes from `pipeline_id` to `pipeline_key = {pipelineName}_v{pipelineVersion}`.
- API contract remains unchanged (`pipelineName`/`pipelineVersion` stay as-is).

Canonical key format: `pipeline_key = {pipelineName}_v{pipelineVersion}`
- all lowercase
- underscore separator before version (example: `array_imputation_v3`)

## Scope
### In scope
- Remove DB definition read paths for:
  - `pipeline_input_definitions`
  - `pipeline_output_definitions`
  - `pipeline_quotas`
- Move pipeline definition fields from `pipelines` to YAML.
- Add `pipeline_runtime_metadata` for runtime-managed fields that are currently in `pipelines`:
  - `pipeline_key`
  - `tool_version`
  - `workspace_billing_project`
  - `workspace_name`
  - `workspace_google_project`
  - `workspace_storage_container_name`
  - `hidden`
- Migrate internal persistence and lookups to `pipeline_key` as the only identifier post-cutover.
- Update `pipeline_runs` to store and use `pipeline_key`; remove `pipeline_id` dependency from runtime read/write paths.
- Perform cutover during planned maintenance window with full in-flight drain before deploy.

### Out of scope
- No API endpoint changes.
- No API request/response schema changes.
- No `common/openapi.yml` changes.
- No generated client regeneration or client migration.

## Key decisions
- **No dual-source mode**: single release cutover from DB-backed definitions to YAML-backed definitions.
- **Maintenance window cutover**: stop accepting new jobs, let old-code jobs complete, then deploy YAML-only readers.
- **No dual-read/dual-write mode**: no long-lived compatibility layer for `pipeline_id` after cutover.
- **Rollback model**: release rollback (not runtime source toggle); DB restore is contingency-only.
- **Behavioral continuity**: keep existing tolerance of mutable runtime metadata between prepare and start.

## Phased plan
### Phase 0 - Guardrails
- Capture ADR constraints: API unchanged, additive migrations only, Stairway compatibility.
- Finalize `pipeline_key` normalization rules.

### Phase 1 - YAML schema and validation
- Extend `pipelines-config.yml` with metadata, inputs, outputs, and quota defaults.
- Add typed config binding and startup validation for uniqueness/completeness.

### Phase 2 - Additive DB migration
- Add `pipeline_key` where needed; backfill from `(pipeline_name, pipeline_version)`.
- Add indexes/uniqueness constraints.
- Create `pipeline_runtime_metadata` with one row per `pipeline_key`.
- Backfill order: `pipeline_key` -> `pipeline_runtime_metadata` -> referential integrity checks.
- Add cutover-ready schema changes for `pipeline_runs` key-based operation and plan removal of `pipeline_id` dependency.

### Phase 3 - Internal refactor and readiness
- Move services to key-centric resolution with YAML-backed definitions.
- Keep API/controller mapping unchanged.
- Persist execution-critical runtime metadata in `pipeline_runtime_metadata` and `pipeline_runs`.
- Remove remaining `pipeline_id` runtime dependencies in services, repositories, and flights.
- Add readiness checks: YAML coverage parity, backfill completeness, historical run resolvability by key.

### Phase 4 - Single-release cutover
- Version flights/steps where working-map contracts change.
- During maintenance window:
  - stop new jobs
  - wait for old in-flight jobs to complete
  - run/verify migrations and backfills
  - verify no incompatible in-flight jobs remain
  - deploy YAML-only + key-only runtime release
  - run immediate health checks

### Phase 5 - Post-cutover validation
- Confirm startup passes with full YAML validation.
- Confirm list/details and prepare/start/result flows work.
- Confirm admin update continues to write `pipeline_runtime_metadata`.
- Confirm runtime reads/writes no longer depend on `pipeline_id`.
- If needed, roll back application release.

### Phase 6 - Cleanup
- Remove dead DB definition read paths.
- Drop obsolete config tables/columns and remaining `pipeline_id` schema after soak period.

## Production cutover runbook checklist
- [ ] **Pre-cutover**: YAML parity confirmed, key backfills complete, runbook approved, on-call staffed.
- [ ] **Window start**: announce maintenance, stop new jobs, drain old in-flight jobs.
- [ ] **Deploy**: run migrations/checks, deploy YAML-only + key-only release, ensure no mixed-source serving.
- [ ] **Validate**: startup, key APIs/flows, admin updates, no runtime reads from definition tables, no `pipeline_id` runtime dependency.
- [ ] **Rollback trigger**: startup failure, endpoint failure, unresolved run/pipeline mappings.
- [ ] **Exit**: stable service, soak started, YAML-only source confirmed, API contract unchanged.

## Testing strategy
- **Unit**: key format/parsing, YAML binding/validation.
- **Migration**: backfill correctness, uniqueness/index constraints.
- **Integration**: prepare/start/result behavior, historical run resolvability by key, no reads from definition tables, no runtime `pipeline_id` usage.
- **Stairway**: old jobs complete pre-cutover; new jobs run post-cutover with key-only working-map fields.
- **Regression/smoke**: unchanged API behavior.

## Risks and mitigations
- In-flight jobs tied to old identifiers -> drain all in-flight jobs before cutover; block deploy until drain is complete.
- Backfill/YAML omissions -> enforce pre-cutover parity, key resolvability checks, and startup fail-fast validation.
- Cutover failure without compatibility fallback -> rehearse release rollback runbook and maintenance-window abort criteria.

## Incremental PR sequence
1. ADR + `pipeline_key` rules + YAML schema scaffold.
2. Additive migrations for `pipeline_key`, indexes, and `pipeline_runtime_metadata`.
3. Internal key-only service/repository/flight refactor + startup validation.
4. Cutover readiness checks and maintenance-window runbook rehearsal.
5. Production cutover release (YAML-only + key-only) + post-cutover validation.
6. Cleanup of legacy DB definition paths/schema, including final `pipeline_id` removal.
