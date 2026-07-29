# scrnaseq-skills-v2

**Status:** v0.1 — Initial scRNA-Seq Pipeline framework. Provisional.

A library of [Claude Code](https://docs.claude.com/claude-code) "skills" for
running single-cell RNA-seq analyses via natural-language prompts. Structured
as a layered hierarchy of primitives, modules, and pipelines that Claude Code
composes when it runs an analysis.

This is a **framework, not a finished product.** It reflects one working
approach to Claude-driven bioinformatics — published in the hope that others
find the structure useful as a starting point for their own skill libraries.

## What this is

If you've used Claude on the web or in an IDE, you've prompted Claude directly
for individual answers. This library is a different pattern: a set of
markdown-based **skills** — reusable, documented analysis units — that
Claude Code loads on demand and composes into full pipelines. You write a
short brief describing what you want; Claude Code reads the relevant skills
and produces a complete analysis.

The library is organized in three layers:

- **`primitives/`** — generic, reusable building blocks (differential
  expression functions, visualization helpers, doublet removal, R environment
  conventions). These are the lowest-level, most stable pieces.
- **`modules/`** — higher-level analyses that compose primitives into
  complete workflows (CellChat, subclustering, pySCENIC regulon inference,
  trajectory analysis, etc.).
- **`pipelines/`** — end-to-end scripts that chain modules together for a
  specific analysis goal (integrate public data, analyze a large dataset).

The `examples/` directory contains worked examples showing how modules are
invoked in practice. The `context/` directory holds cross-cutting reference
material (color palettes, known public atlases, lab conventions).

## Who this is for

Researchers comfortable with scRNA-seq analysis (Seurat, standard workflows)
who want to try an agent-driven pipeline approach. No prior Claude Code
experience needed — `GETTING_STARTED.md` walks through setup.

## What state this is in

**v0.1 — Initial framework.** This is a first public release of a working
internal library. Expect:

- Solid, well-tested primitives and modules for common workflows
- Rough edges in less-used skills
- Documentation gaps in places
- Ongoing evolution — this is not a stable API

See `STATUS.md` for a per-directory maturity map.

## Prerequisites

- Claude Code installed and configured — see [Anthropic's docs](https://docs.claude.com/claude-code)
- A conda environment named `r-env` with R 4.5+, Seurat v5, and standard scRNA-seq packages (setup notes in `GETTING_STARTED.md`)
- For pySCENIC regulon inference: a separate conda environment named `scenicenv`
- Your own scRNA-seq data (10x cellranger output is the assumed format)

## Getting started

1. Clone this repo:
   ```
   git clone https://github.com/Ryannachmand/scrnaseq-skills-v2.git ~/claude-skills-v2
   ```
   The path `~/claude-skills-v2` is what several files reference — clone
   there or update the paths accordingly.
2. Read `GETTING_STARTED.md`.
3. Work through the example prompt at the end of that file, against your
   own data.

## Layout

```
claude-skills-v2/
├── SKILL.md                  ← entry point Claude reads first
├── CONVENTIONS.md            ← naming, file paths, config style
├── PROJECT_CLAUDE_TEMPLATE.md ← template for per-project CLAUDE.md files
│
├── primitives/               ← generic building blocks
├── modules/                  ← composed analyses
├── pipelines/                ← end-to-end workflows
├── examples/                 ← worked examples
├── context/                  ← reference material (palettes, atlases)
│
├── README.md                 ← this file
├── GETTING_STARTED.md        ← onboarding walkthrough
└── STATUS.md                 ← per-directory maturity
```

## Feedback

This is a v0.1 public release; the structure will evolve. Issues and PRs
welcome — no expectation of stable behavior yet.
