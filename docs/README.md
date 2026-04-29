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

## Manuel updated of LLM documentation files

For now we need to get the API documentation in markdown

> We need to update these files manually on PR-branches.
> On the main branch these files are updated automatically by the CI for a commit.

```
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
# build docs
sphinx-build -n -W --keep-going -b html ./ ./_build/
```
