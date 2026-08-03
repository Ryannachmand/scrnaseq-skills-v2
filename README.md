# scrnaseq-skills-v2

**Status:** v0.1 — Initial scRNA-seq framework. Experimental and under active development.

A library of **Claude Code** skills (natural-language instructions and reusable R code) for performing single-cell RNA-seq analyses through natural-language prompts.

This repository organizes computational biology expertise into modular, reusable analysis units that Claude Code can compose into complete workflows. The goal is not to automate scientific judgment, but to provide a structured, reusable foundation for AI-assisted computational biology. 

This repository reflects one working approach to organizing AI-driven bioinformatics workflows. It is an experimental framework rather than a finished product, and is intended both as a practical tool and as a reference for researchers interested in building similar AI-assisted analysis systems.

---

## Why this exists

Modern single-cell RNA-seq analyses combine many interconnected computational steps: quality control, normalization, integration, clustering, annotation, downstream biological analyses, and publication-quality visualization. While AI coding agents can accelerate these workflows, they perform best when supplied with explicit, reusable domain knowledge rather than relying entirely on ad hoc prompting.

This repository explores one approach to that problem by encoding analysis workflows as modular "skills" that can be reused across projects. Instead of rewriting prompts or copying analysis code between repositories, computational workflows are organized into reusable building blocks that an AI coding agent can assemble as needed.

---

## What this is

If you've used Claude through the web interface or an IDE, you've probably interacted with it through individual prompts. This repository supports a different workflow.

It contains markdown-based **skills**—reusable, documented analysis units—that Claude Code loads on demand while planning an analysis. Some skills contain reusable R code templates, while others encode workflow guidance, conventions, or biological analysis strategies. Given a short description of the desired analysis, Claude Code selects the relevant skills and composes them into a complete workflow.

The library is organized into three layers:

* **`primitives/`** — low-level reusable building blocks such as differential expression, visualization utilities, quality control, doublet detection, and common R conventions.
* **`modules/`** — higher-level biological analyses that combine primitives into complete workflows (CellChat, pySCENIC, trajectory analysis, subclustering, etc.).
* **`pipelines/`** — end-to-end analyses that combine multiple modules for common research tasks such as integrating public datasets or performing comprehensive downstream analyses.

The `examples/` directory contains worked examples demonstrating how these components are used together. The `context/` directory contains shared reference material including color palettes, public atlases, and laboratory conventions.

---

## How it works

A typical workflow is:

1. Describe the biological question or desired analysis in natural language.
2. Claude Code loads the relevant skills from this repository.
3. Appropriate primitives, modules, and pipelines are selected.
4. Claude generates and executes customized analysis code for the specific dataset.
5. Figures, tables, and analysis outputs are produced for downstream interpretation.

The researcher remains responsible for evaluating biological conclusions; this framework is intended to standardize computational workflow execution rather than replace scientific decision-making.

---

## Who this is for

Researchers who:

* routinely analyze single-cell RNA-seq data using Seurat or related tools
* want to incorporate AI coding agents into their computational workflows
* value reusable and reproducible analysis infrastructure
* are interested in experimenting with modular AI-assisted bioinformatics

No prior Claude Code experience is required. See `GETTING_STARTED.md` for installation and setup.

---

## Current status

**v0.1** is the first public release of a working internal framework.

Current strengths include:

* reusable analysis primitives for common workflows
* modular organization of biological analyses
* support for many standard downstream scRNA-seq analyses
* documented conventions for AI-assisted workflow execution

Current limitations include:

* uneven maturity across advanced analyses
* documentation gaps in some modules
* evolving interfaces and conventions
* no guarantee of API stability

See `STATUS.md` for a detailed overview of repository maturity.

---

## Prerequisites

* Claude Code installed and configured
* A conda environment named `r-env` containing R 4.5+, Seurat v5, and standard scRNA-seq packages (see `GETTING_STARTED.md`)
* A separate `scenicenv` environment for pySCENIC analyses
* Your own scRNA-seq dataset (10x Cell Ranger output is the assumed input format)

---

## Getting started

```bash
git clone https://github.com/Ryannachmand/scrnaseq-skills-v2.git ~/claude-skills-v2
```

Then:

1. Read `GETTING_STARTED.md`
2. Configure the required environments
3. Work through the example project using your own dataset

---

## Repository layout

```text
scrnaseq-skills-v2/
├── SKILL.md
├── CONVENTIONS.md
├── PROJECT_CLAUDE_TEMPLATE.md
│
├── primitives/
├── modules/
├── pipelines/
├── examples/
├── context/
│
├── README.md
├── GETTING_STARTED.md
└── STATUS.md
```

---

## Feedback

This is an experimental framework and will continue to evolve as new workflows are developed and existing skills mature.

Issues, suggestions, and pull requests are welcome.
