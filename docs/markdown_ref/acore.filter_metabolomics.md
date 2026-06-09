# acore.filter_metabolomics package

Module for filtering metabolomics feature table.

### filter_by_missingness(data: DataFrame, percent: [int](https://docs.python.org/3/library/functions.html#int) = 80, method: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'classic', samples: [list](https://docs.python.org/3/library/stdtypes.html#list) | [None](https://docs.python.org/3/library/constants.html#None) = None, groups: [dict](https://docs.python.org/3/library/stdtypes.html#dict) | [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None) = None)

Implementation of the 80%-rule.

If there are more than 20% of values (intensities) missing for one feature,
this feature will get removed.

* **Parameters:**
  * **data** – pandas data frame with samples as rows and features as columns.
  * **percent** – percentage chosen for filtering. The default is 80%, meaning that
    at least 80% of the values of every feature need to be present in order for this
    feature to be retained.
  * **method** – str that is either “classic” or “modified”.
    If “classic”, all samples are considered for each feature. Samples are taken from the
    “samples” parameter and should not include controls or QCs.
    If “modified”, conditions are separated when calculating the percentage of
    missingness. A feature is retained if at least 

    ```
    ``
    ```

    percent\`\`% of its values are
    present in ANY one condition. This allows condition-specific features (e.g.
    present in treatment but missing in control) to be retained.
  * **samples** – list of row index labels (from data.index) identifying the biological
    sample rows, e.g. [“S1”, “S2”, “S3”]. Required when method=”classic”. Should not
    include control or QC samples.
  * **groups** – 

    required when method=”modified”, ignored otherwise. Can be either:
    - A dict mapping condition name to a list of row index labels belonging to that
      condition, e.g. `{"treatment": ["S1", "S2", "S3"], "control": ["S4", "S5"]}`.
      QCs and blanks are excluded by simply not including them in the dict.
    - A str naming a column in `data` whose values define the condition for each row,
      e.g. `"sample collection"` if rows carry values like `"Berlin"`, `"Copenhagen"`, `"London"`.
      Every unique value in that column becomes a condition group containing all rows
      with that value. When using this option, make sure to not include any other metadata
      columns in the data frame.

### filter_cv(data: DataFrame, samples: [list](https://docs.python.org/3/library/stdtypes.html#list), qcs: [list](https://docs.python.org/3/library/stdtypes.html#list))

Implementation of coefficient of variation (CV)-based filtering.

Features are removed when their CV across biological samples is smaller than their CV
across QC samples, meaning analytical noise exceeds biological variability.

* **Parameters:**
  * **data** – pandas data frame with samples as rows and features as columns.
  * **samples** – list of row index labels (from data.index) identifying the
    biological sample rows, e.g. [“S1”, “S2”, “S3”].
  * **qcs** – list of row index labels identifying the quality control rows,
    e.g. [“QC1”, “QC2”, “QC3”].

### filter_blanks(data: DataFrame, blanks: [list](https://docs.python.org/3/library/stdtypes.html#list), samples: [list](https://docs.python.org/3/library/stdtypes.html#list), threshold: [float](https://docs.python.org/3/library/functions.html#float) = 0.5)

Filtering out features that show up in the blanks control.

The mean intensity scores are calculated per-feature within the
blanks and the samples. If the ratio of a feature’s mean intensity in the blanks
to its mean intensity in the samples is more than half (per default), the feature
gets removed. It is assumed to have potentially contaminated the instrument, so
the measurements in the samples cannot be trusted to be biologically relevant.

* **Parameters:**
  * **data** – pandas DataFrame containing data with samples as rows and features as columns
  * **blanks** – list of row index labels (from data.index) identifying the blanks
    measurement rows, e.g. [“Blank1”, “Blank2”]
  * **samples** – list of row index labels (from data.index) identifying the biological
    sample rows, e.g. [“S1”, “S2”, “S3”]
  * **threshold** – optional ratio used as a threshold to determine whether the detected
    intensities in blanks are too high in comparison with sample intensities.
    Defaults to 0.5, but can be adjusted based on data and stringency.
