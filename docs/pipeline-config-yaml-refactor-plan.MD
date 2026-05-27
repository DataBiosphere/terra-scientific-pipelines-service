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
- Keep `pipeline_runs` with both `pipeline_id` and `pipeline_key` for compatibility.
- Preserve in-flight Stairway jobs via versioned flights/steps.

### Out of scope
- No API endpoint changes.
- No API request/response schema changes.
- No `common/openapi.yml` changes.
- No generated client regeneration or client migration.

## Key decisions
- **No dual-source mode**: single release cutover from DB-backed definitions to YAML-backed definitions.
- **Maintenance window cutover**: stop accepting new jobs, let old-code jobs complete, then deploy YAML-only readers.
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

### Phase 3 - Internal refactor and readiness
- Move services to key-centric resolution with YAML-backed definitions.
- Keep API/controller mapping unchanged.
- Persist execution-critical runtime metadata in `pipeline_runtime_metadata` and `pipeline_runs`.
- Add readiness checks: YAML coverage parity, backfill completeness, historical run resolvability.

### Phase 4 - Single-release cutover
- Version flights/steps where working-map contracts change.
- During maintenance window:
  - stop new jobs
  - wait for old in-flight jobs to complete
  - run/verify migrations and backfills
  - verify no incompatible in-flight jobs remain
  - deploy YAML-only definition reader release
  - run immediate health checks

### Phase 5 - Post-cutover validation
- Confirm startup passes with full YAML validation.
- Confirm list/details and prepare/start/result flows work.
- Confirm admin update continues to write `pipeline_runtime_metadata`.
- If needed, roll back application release.

### Phase 6 - Cleanup
- Remove dead DB definition read paths.
- Drop obsolete config tables/columns after soak period.

## Production cutover runbook checklist
- [ ] **Pre-cutover**: YAML parity confirmed, backfills complete, runbook approved, on-call staffed.
- [ ] **Window start**: announce maintenance, stop new jobs, drain old in-flight jobs.
- [ ] **Deploy**: run migrations/checks, deploy YAML-only release, ensure no mixed-source serving.
- [ ] **Validate**: startup, key APIs/flows, admin updates, no runtime reads from definition tables.
- [ ] **Rollback trigger**: startup failure, endpoint failure, unresolved run/pipeline mappings.
- [ ] **Exit**: stable service, soak started, YAML-only source confirmed, API contract unchanged.

## Testing strategy
- **Unit**: key format/parsing, YAML binding/validation.
- **Migration**: backfill correctness, uniqueness/index constraints.
- **Integration**: prepare/start/result behavior, historical run resolvability, no reads from definition tables.
- **Stairway**: old jobs complete pre-cutover; new jobs run post-cutover.
- **Regression/smoke**: unchanged API behavior.

## Risks and mitigations
- In-flight jobs tied to old identifiers -> version flights/steps and drain during maintenance window.
- Backfill/YAML omissions -> enforce pre-cutover parity and startup fail-fast validation.
- Cutover failure without dual-source fallback -> rehearse release rollback runbook.

## Incremental PR sequence
1. ADR + `pipeline_key` rules + YAML schema scaffold.
2. Additive migrations for `pipeline_key`, indexes, and `pipeline_runtime_metadata`.
3. Internal key-centric service refactor + startup validation.
4. Stairway compatibility updates.
5. YAML source-of-truth implementation + cutover prep.
6. Production cutover release + post-cutover validation.
7. Cleanup of legacy DB definition paths/schema.
