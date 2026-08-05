<!--
Please complete the following sections when you submit your pull request.
Note that text within html comment tags will not be rendered.

Have a look at the CONTRIBUTING.md file in the root of the repository. -->

### Summary

<!-- Describe the problem you're trying to fix in this pull request. Please reference any related issue and use fixes/close to automatically close them, if pertinent. For example: "Fixes #58", or "Addresses (but does not close) #238". -->

Fixes #<NUM>

...

### List of changes proposed in this PR (pull-request)

<!-- We suggest using bullets (indicated by * or -) and filled checkboxes [x] here -->

- [ ] ...
- [ ] ...

### Checks

- [ ] Look at [CONTRIBUTING.md](https://github.com/Multiomics-Analytics-Group/acore/blob/HEAD/CONTRIBUTING.md) for guidance on how to contribute.

In case you add a new module or update one, please check the following tasks:

- [ ] The folder names defines the name of the subpackage (module)
- [ ] Add the user-facing functions to the `__init__.py` in the new folder, so that
      they are available when the subpackage is imported.
- [ ] Create **Pandera schema** in a file with subpackage name in the `src/acore/types` folder.
      Optimal is to have only one output schema of results per subpackage or module.
- [ ] Add a relatively small and public dataset to the `data` folder, or reuse an existing one for testing
- [ ] Create an **api example jupyter notebook** in the `docs/api_examples_module` folder with that data
- [ ] Use **jupytext** to sync the Jupyter notebook with a Python script
- [ ] Update `index.md` file in the `docs` folder with the new example
- [ ] Create **test script** in the `/tests` folder with the name of the subpackage or module
      using pytest or unittests to test your new functionality.
