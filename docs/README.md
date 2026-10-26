# Docs creation

In order to build the docs you need to

1. install sphinx and additional support packages
2. build the package reference files
3. run sphinx to create a local html version

The documentation is build using readthedocs automatically.

Install the docs dependencies of the package (as speciefied in toml):

```bash
# in main folder
pip install .[docs]
```

## Updating LLM documentation files

The API reference in `docs/markdown_ref` is updated automatically by the CI workflow.
If you need to regenerate these files locally, or if the automation is unavailable for your branch,
run the following commands manually:

```bash
# we are only interested in the reference files
sphinx-apidoc --force --implicit-namespaces --module-first -o reference ../
sphinx-build -n -W --keep-going -b markdown ./ ./_build_markdown
mv _build_markdown/reference ./markdown_ref
```

## Build docs using Sphinx command line tools

Command to be run from `path/to/docs`, i.e. from within the `docs` package folder:

Options:

- `--separate` to build separate pages for each (sub-)module

```bash
# pwd: docs
# apidoc
sphinx-apidoc --force --implicit-namespaces --module-first -o reference ../src/acore
# download other resources
python fetch_files.py
# parse main repository README
python split_readme.py
# build docs (-j auto executes the notebooks in parallel)
sphinx-build -n -W --keep-going -j auto -b html ./ ./_build/
```

## Notebooks depending on external APIs

Notebooks which mainly download data from external APIs (PRIDE, journal websites) are
not executed during the docs build (see `nb_execution_excludepatterns` in `conf.py`).
They are executed in the separate `Test API notebooks` GitHub Actions workflow
(`.github/workflows/api_notebooks.yml`) together with the API tests in `tests_api`,
weekly and on changes to `src/acore/io`.

The enrichment analysis example loads cached UniProt annotations from
`example_data/PXD010372/processed/annotations_100.csv` and only queries UniProt if
that file is not available. In that case the fetched annotations are saved to this
path (and uploaded as an artifact by the `Test API notebooks` workflow).
