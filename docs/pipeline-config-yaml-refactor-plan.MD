# Pipeline Config YAML Refactor Plan [ARCHIVED - SEE COMPLETION SUMMARY]

⚠️ **This is the original detailed refactor plan document. For current status and remaining work, see [`pipeline-config-yaml-refactor-completion.MD`](./pipeline-config-yaml-refactor-completion.MD).**

## Executive Summary

**Refactor Status: ~95% Complete** ✅

The pipeline definitions refactor from DB tables to YAML-backed configuration with `pipeline_key` identifiers is mostly complete:
- ✅ Phase 1: YAML schema with full definitions and comprehensive validation
- ✅ Phase 2: Database schema additions (pipeline_runtime_metadata, pipeline_key)
- ✅ Phase 3: Service layer refactor (PipelinesService YAML-backed, pipeline_key-centric)
- 📝 Remaining: Definition-centric test refactor + legacy cleanup

**No maintenance window was needed.** We used additive migrations + hybrid identifier approach (both `pipeline_id` and `pipeline_key`) to deploy without downtime.

---

## What's Implemented

### Canonical Identifier
- **`pipeline_key`**: `{pipelineName}_v{pipelineVersion}` (lowercase, e.g., `array_imputation_v2`)
- Used throughout service layer, Stairway flights, and database
- Replaces `pipeline_id` for all new operations

### Configuration & Validation
- `pipelines-config.yml` fully populated with all pipeline definitions
- `PipelineConfigurations` provides typed Spring Boot binding
- `PipelineConfigurationsTest` validates YAML (80+ test cases covering structure, inputs, outputs, quotas)

### Database
- `pipeline_runtime_metadata` table created (PK: pipeline_key)
- `pipeline_runs` enhanced: `pipeline_key` added (non-nullable), `pipeline_id` remains (nullable)
- All migrations applied; backfills completed

### Service Layer
- `PipelinesService` reads definitions from YAML + merges with `PipelineRuntimeMetadata`
- API controllers unchanged (backward compatible)
- No breaking changes to API contract

---

## What's Remaining

For details on remaining work and migration checklists, see **[`pipeline-config-yaml-refactor-completion.MD`](./pipeline-config-yaml-refactor-completion.MD)**.

Quick summary:
- 📝 Refactor decision-test fixture from `Pipeline.java` entities to YAML fixtures
- 🧹 Remove legacy DB definition table read paths post-validation
- ♻️ Drop obsolete tables after soak period

---

## Key Implementation Decisions That Held

| Decision | Rationale | Status |
|----------|-----------|--------|
| **API contract unchanged** | Minimize client impact; hide refactor complexity | ✅ Held |
| **Additive migrations** | No downtime; backward-compatible schema evolution | ✅ Held |
| **Dual identifiers (transitional)** | `pipeline_id` + `pipeline_key` simplifies migration | ✅ Held |
| **YAML as definition source of truth** | Configuration-as-code; easier version control + reviews | ✅ Held |
| **Comprehensive validation** | Prevent misconfigured YAML reaching production | ✅ Held |
| **No compatibility layer** | Direct pipeline_key resolution; reduces code complexity | ✅ Held |

## What Was Different from Original Plan

Original plan assumed a maintenance-window cutover with drain-and-redeploy. In practice:
- We deployed with **hybrid identifiers** (both `pipeline_id` and `pipeline_key`) to avoid downtime
- **Validation-first** approach: deployed comprehensive tests before production code
- **Parallel execution** of phases where possible to reduce feedback cycles
- Avoided the maintenance window entirely through careful migration design

## Key Resources

| Resource | Purpose | Location |
|----------|---------|----------|
| **YAML Config** | Pipeline definitions source of truth | `service/src/main/resources/pipelines-config.yml` |
| **Config Classes** | Spring Boot binding; typed representation | `service/src/main/java/bio/terra/pipelines/app/configuration/internal/PipelineConfigurations.java` |
| **Validation Test** | Comprehensive YAML validation (80+ tests) | `service/src/test/java/bio/terra/pipelines/configuration/internal/PipelineConfigurationsTest.java` |
| **Service** | YAML-backed, pipeline_key-centric resolution | `service/src/main/java/bio/terra/pipelines/service/PipelinesService.java` |
| **DB Entities** | PipelineRuntimeMetadata, PipelineRun | `service/src/main/java/bio/terra/pipelines/db/entities/` |
| **Completion Doc** | Current status & remaining work | `docs/pipeline-config-yaml-refactor-completion.MD` |

---

**For detailed remaining work items, migration checklist, and medium/long-term cleanup strategy, see [`pipeline-config-yaml-refactor-completion.MD`](./pipeline-config-yaml-refactor-completion.MD).**


