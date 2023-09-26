========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |github-actions|
        | |codecov|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|
.. |docs| image:: https://readthedocs.org/projects/acore/badge/?style=flat
    :target: https://acore.readthedocs.io/
    :alt: Documentation Status

.. |github-actions| image:: https://github.com/albsantosdel/acore/actions/workflows/github-actions.yml/badge.svg
    :alt: GitHub Actions Build Status
    :target: https://github.com/albsantosdel/acore/actions

.. |codecov| image:: https://codecov.io/gh/albsantosdel/acore/branch/main/graphs/badge.svg?branch=main
    :alt: Coverage Status
    :target: https://app.codecov.io/github/albsantosdel/acore

.. |version| image:: https://img.shields.io/pypi/v/acore.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/acore

.. |wheel| image:: https://img.shields.io/pypi/wheel/acore.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/acore

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/acore.svg
    :alt: Supported versions
    :target: https://pypi.org/project/acore

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/acore.svg
    :alt: Supported implementations
    :target: https://pypi.org/project/acore

.. |commits-since| image:: https://img.shields.io/github/commits-since/albsantosdel/acore/v0.1.0.svg
    :alt: Commits since latest release
    :target: https://github.com/albsantosdel/acore/compare/v0.1.0...main



.. end-badges

A Python package with statistical functions to analyse multimodal molecular data

* Free software: GNU Lesser General Public License v3 or later (LGPLv3+)

Installation
============

::

    pip install acore

You can also install the in-development version with::

    pip install https://github.com/albsantosdel/acore/archive/main.zip


Documentation
=============


https://acore.readthedocs.io/


Development
===========

To run all the tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
