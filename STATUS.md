# STATUS

Per-directory maturity map for `scrnaseq-skills-v2`, as of the initial
public release (v0.1).

**Overall status: v0.1 — Initial framework.** Structure is stable enough to
build on but expect ongoing evolution. This document exists to be honest
about which parts of the library you can lean on and which parts will
change.

---

## Legend

- **Stable** — well-tested, unlikely to change shape. Safe to depend on.
- **Working** — functional and used regularly, but may have rough edges,
  documentation gaps, or minor evolution over time.
- **Provisional** — usable but expect changes, gaps, or friction. Read the
  code before depending on the behavior.

---

## By directory

### `primitives/` — Stable

Low-level reusable functions (differential expression, visualization,
doublet removal, R environment conventions). Most heavily audited layer of
the library. Function signatures and behavior are unlikely to change
without a version bump.

### `modules/` — Working

Higher-level analyses composed from primitives (CellChat, subclustering,
trajectory, pySCENIC, metabolic profiling, cross-dataset dotplots, etc.).
Each module has been validated on real datasets. Expect:

- Occasional edge cases the current implementation doesn't handle
- Module-internal parameter defaults that may not suit every dataset
- Some modules more polished than others; `cellchat.md` and
  `de_comprehensive_csv.md` are among the most-used, others less so

Read the module file before running to understand what parameters matter
for your data.

### `pipelines/` — Working

Two end-to-end pipelines are shipped:

- `pipelines/large_dataset/` — analyze one in-house dataset end-to-end
  (loading, QC, integration, clustering, annotation, optional subset
  analysis, optional downstream modules). Validated on typical multi-sample
  10x runs.
- `pipelines/integrate_public_data/` — integrate in-house data with a
  public reference atlas (e.g., Tabula Sapiens). Validated on the atlases
  registered in `context/known_atlases.md`.

Pipelines are the primary user-facing entry points. Brief templates
(`brief_template.txt`) document all required and optional fields.

### `examples/` — Provisional

Worked examples showing how each module is invoked in practice. Files
follow the `<module>_example_A.md` naming convention.

- **Purpose is pedagogical.** Examples are meant to be read to understand
  invocation patterns, not run verbatim.
- **Values marked `<PLACEHOLDER>` or `# REPLACE:` must be customized**
  for your data before use.
- **Coverage is not exhaustive** — not every module has a worked example
  yet. Expect additions over time.

The `tabula_sapiens_*.md` examples reference the public Tabula Sapiens
atlas and can be adapted directly if you use that atlas as a reference.

### `context/` — Working

Cross-cutting reference material loaded by pipelines and modules:

- `color_palettes.md` — named palettes referenced by group_col and
  subtype_col in briefs
- `known_atlases.md` — registry of recognized public reference atlases
  with canonical parameters
- `lab_context.md` — conventions for column naming, cross-dataset
  ordering rules, canonical marker gene sets, and analysis defaults
- `validated_examples.yaml` — machine-readable index of validated example
  runs (partial coverage; entries added as new pipelines are run)

Content evolves as new datasets and analyses are validated. Structure is
stable.

---

## Root-level files

- `SKILL.md` — entry point Claude reads first. **Stable.**
- `CONVENTIONS.md` — naming, path, and config conventions. **Stable.**
- `PROJECT_CLAUDE_TEMPLATE.md` — template for per-project `CLAUDE.md`
  files. **Stable**, but references the path `~/claude-skills-v2/`
  explicitly — update if you clone elsewhere.

---

## Roadmap

No formal roadmap for v0.1. Priorities for future versions include broader
example coverage, additional pipeline entry points for common downstream
questions, and refinements to module parameter defaults based on
feedback.

Issues and PRs welcome on GitHub.
