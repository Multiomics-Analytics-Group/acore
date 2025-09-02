## ![ACore Logo](https://raw.githubusercontent.com/Multiomics-Analytics-Group/acore/HEAD/docs/images/logo/acore_logo.svg)

<p align="center">
   ACore is Python package with statistical functions to analyse multimodal molecular data
</p>

| Information           | Links                                                                                                                                                                                                                                                                                                                                                       |
| :-------------------- | :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Package**           | [![PyPI Latest Release](https://img.shields.io/pypi/v/acore.svg)][acore-pypi] [![Supported versions](https://img.shields.io/pypi/pyversions/acore.svg)][acore-pypi] [![Supported implementations](https://img.shields.io/pypi/implementation/acore.svg)][acore-pypi] [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)][gpl-license] |
| **Documentation**     | [![View - Documentation](https://img.shields.io/badge/view-Documentation-blue?style=flat)][acore-docs] [![made-with-sphinx-doc](https://img.shields.io/badge/Made%20with-Sphinx-1f425f.svg)](https://www.sphinx-doc.org/) ![Docs](https://readthedocs.org/projects/analytics-core/badge/?style=flat) [![CC BY 4.0][cc-by-shield]][acore-license]            |
| **Build**             | [![CI](https://github.com/Multiomics-Analytics-Group/acore/actions/workflows/tox-gha.yml/badge.svg)][ci-gh-action]                                                                                                                                                                                                                                          |
| **Discuss on GitHub** | [![GitHub issues](https://img.shields.io/github/issues/Multiomics-Analytics-Group/acore)][issues] [![GitHub pull requests](https://img.shields.io/github/issues-pr/Multiomics-Analytics-Group/acore)][pulls]                                                                                                                                                |

## Table of contents:

- [About the project](#about-the-project)
- [Installation](#installation)
- [Documentation](#documentation)
- [License](#license)
- [Contributing](#contributing)
- [Credits and acknowledgements](#credits-and-acknowledgements)
- [Contact and feedback](#contact-and-feedback)

## About the project

ACore is Python package with statistical functions to analyse multimodal molecular data. It is part of a broader ecosystem of tools for multi-omics analysis, working in conjunction with [VueCore][acore] and [VueGen][vuegen] to enable end-to-end data processing, visualization, and reporting.

## Installation

> [!TIP]
> It is recommended to install ACore inside a virtual environment to manage depenendencies and avoid conflicts with existing packages. You can use the virtual environment manager of your choice, such as `poetry`, `conda`, or `pipenv`.

### Pip

ACore is available on [PyPI][acore-pypi] and can be installed using pip:

```bash
pip install acore
```

You can also install the package for development by cloning this repository and running the following command:

> [!WARNING]
> We assume you are in the root directory of the cloned repository when running this command. Otherwise, you need to specify the path to the `acore` directory.

```bash
pip install -e '.[dev]'
```

To run all tests, you can use the following command:

```bash
pytest
```

## Documentation

ACore's documentation is hosted on [Read the Docs][acore-docs]. It includes installation instructions, various tutorials, configuration options, and the API reference.

## License

The code in this repository is licensed under the **GNU General Public License v3**, allowing you to use, modify, and distribute it freely as long as you include the original copyright and license notice.

The documentation and other creative content are licensed under the **Creative Commons Attribution 4.0 International (CC BY 4.0) License**, meaning you are free to share and adapt it with proper attribution.

Full details for both licenses can be found in the [LICENSE][acore-license] file.

## Contributing

ACore is an open-source project, and we welcome contributions of all kinds via GitHub issues and pull requests. You can report bugs, suggest improvements, propose new features, or implement changes. Please follow the guidelines in the [CONTRIBUTING](CONTRIBUTING.rst) file to ensure that your contribution is easily integrated into the project.

## Credits and acknowledgements

- ACore was developed by the [Multiomics Network Analytics Group (MoNA)][Mona] at the [Novo Nordisk Foundation Center for Biosustainability (DTU Biosustain)][Biosustain].

## Contact and feedback

We appreciate your feedback! If you have any comments, suggestions, or run into issues while using ACore, feel free to [open an issue][new-issue] in this repository. Your input helps us make ACore better for everyone.

[acore-pypi]: https://pypi.org/project/acore/
[acore-license]: https://github.com/Multiomics-Analytics-Group/acore/blob/main/LICENSE.md
[acore-docs]: https://analytics-core.readthedocs.io/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
[gpl-license]: https://www.gnu.org/licenses/gpl-3.0
[ci-gh-action]: https://github.com/Multiomics-Analytics-Group/acore/actions/workflows/tox-gha.yml
[issues]: https://github.com/Multiomics-Analytics-Group/acore/issues
[pulls]: https://github.com/Multiomics-Analytics-Group/acore/pulls
[vuegen]: https://github.com/Multiomics-Analytics-Group/vuegen
[vuecore]: https://github.com/Multiomics-Analytics-Group/vuecore
[Mona]: https://multiomics-analytics-group.github.io/
[Biosustain]: https://www.biosustain.dtu.dk/
[new-issue]: https://github.com/Multiomics-Analytics-Group/acore/issues/new
