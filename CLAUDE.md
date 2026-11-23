# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in
this repository.

## Project

`acore` is a Python package (src layout, `src/acore`) of statistical functions for
multi-omics data: normalization, imputation, differential regulation, enrichment, batch
correction, decomposition, etc. It is the analysis part of an ecosystem with `vuecore`
(visualization) and `vuegen` (reports). Version comes from git tags via `setuptools_scm`.
Python >=3.10.

## Commands

```bash
pip install -e '.[dev]'          # dev install (add '.[all]' for optional network/KM/TDA modules, '.[docs]' for docs)
pytest tests                      # CI runs `python -m pytest tests`
pytest tests/test_normalization.py::test_name   # single test
ruff check src                    # lint (CI)
black . && isort .                # format (CI checks black + isort, profile=black)
pytest tests_api                  # tests hitting external APIs (UniProt etc.); run separately, weekly CI only
```

Pre-commit (`.pre-commit-config.yaml`) runs black, isort and `ruff check` on `src/`,
mirroring CI; set up once per clone with `pip install pre-commit && pre-commit install`,
run manually with `pre-commit run --all-files`. Don't bypass it with `--no-verify`.

Docs: `pip install .[docs]`, then `sphinx-build` from `docs/` (see `docs/README.md`).
Notebooks in `docs/api/` are paired with `.py` (jupytext percent format, configured in
`docs/api/jupytext.toml`) — edit/commit both. Notebooks that depend on external APIs are
excluded from docs execution (`nb_execution_excludepatterns` in `docs/conf.py`) and
tested by `.github/workflows/api_notebooks.yml`.

## Architecture

- Each analysis area is a subpackage of `src/acore/` (e.g. `normalization`,
  `imputation_analysis`, `differential_regulation`, `enrichment_analysis`,
  `multiple_testing`, `batch_correction`, `decomposition`, `exploratory_analysis`). Only
  the core ones are imported in `acore/__init__.py`; modules needing optional
  dependencies (`network_analysis`, `kaplan_meier_analysis`, `tda_analysis`,
  `publications_analysis`, ...) are not, and rely on the `all` extra.
- `acore/types/` holds `pandera` `DataFrameModel` schemas for the result DataFrames of
  analyses (e.g. `EnrichmentAnalysisSchema` in `types/enrichment_analysis.py`). When
  adding or changing result columns in an analysis (e.g. odds ratio columns in
  enrichment), update the matching schema too.
- `enrichment_analysis/`: `annotate.py` and `permutations_enrichment.py` plus
  `statistical_tests/` (`fisher.py` for over-representation, `kolmogorov_smirnov.py` for
  ranked enrichment). The current branch work is on over-representation odds ratios
  (`odds_ratio`, `odds_ratio_lower`; lower bound is for ranking/plotting only).
- Dependency pins carry ToDo notes in `pyproject.toml` (`pandas<3`, `pingouin<0.6`,
  `scipy>=1.10` for the Fisher odds-ratio confidence interval); check them before
  bumping.
- `skills/acore/` contains an agent skill (SKILL.md + references) describing the acore
  API; keep it consistent when public API changes. `docs/skills.md` documents it.
- `docs/` is Sphinx (myst-nb); `docs/reference` is API reference, `docs/api` and
  `docs/analyses` hold executable example notebooks, `docs/sections_readme` and
  `split_readme.py` derive docs pages from the README.

## Conventions

- Commit messages in this repo use gitmoji prefixes (e.g. `:art: format`, `:wrench:`,
  `🚧`).
- PRs target `main`; a PR template is in `.github/PULL_REQUEST_TEMPLATE.md`.
