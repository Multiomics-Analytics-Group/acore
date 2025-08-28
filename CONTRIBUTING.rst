.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/Multiomics-Analytics-Group/acore/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

acore could always use more documentation, whether as part of the
official acore docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/Multiomics-Analytics-Group/acore/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `acore` for local development. Consider an 
advanced editor to help you with some of the common steps described below, e.g. 
`VSCode <https://code.visualstudio.com/docs/introvideos/basics>`_.

1. Fork the `acore` repo on GitHub.
2. Clone your fork locally::

    $ git clone https://github.com/Multiomics-Analytics-Group/acore.git

3. Install your local copy into a virtual environment. Assuming you have Python available 
   on your system, this can be done using `venv`. Alternatives are conda environments
   or uv to create and manage virtual environments::

    $ cd acore/
    $ python -m venv .env
    $ source .env/bin/activate
    $ pip install -e .[dev]

If you work on a Windows shell, see the docs for instructions: 
`How venvs work <https://docs.python.org/3/library/venv.html#how-venvs-work>`_

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the
   tests, including testing other Python versions with tox::

    $ flake8 acore tests
    $ python setup.py test or pytest
    $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should pass the workflows on GitHub.

See for example the PR-Template for a module: 
`Add module PR template <https://github.com/Multiomics-Analytics-Group/acore/blob/main/.github/PULL_REQUEST_TEMPLATE/module.md>`_

Tips
----

To run a subset of tests::

$ pytest tests


Deploying
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in HISTORY.rst).
Then run create a new `GitHub release <https://github.com/Multiomics-Analytics-Group/acore/releases>`_.
The code will then be deployed to PyPI if the tests pass.
