<!--
Please complete the following sections when you submit your pull request.
Note that text within html comment tags will not be rendered.
-->

### Description

Add a subpackage to core.

### Tasks Checklist

- [ ] The folder names defines the name of the subpackage (module)
- [ ] Add the user-facing functions to the `__init__.py` in the new folder, so that
      they are available when the subpackage is imported.
- [ ] Create **Pandera schema** in a file with subpackage name in the `src/acore/types` folder.
      Optimal is to have only one output schema of results per subpackage or module.
- [ ] Create an **api example jupyter notebook** in the `docs/api_examples folder`
- [ ] Use **jupytext** to sync the Jupyter notebook with a Python script
- [ ] Update `index.md` file in the `docs` folder with the new example
- [ ] Create **test script** in the `/test` folder with the name of the subpackage or module
      using pytest or unittests to test your new functionality.
