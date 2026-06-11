# acore.drift_correction package

Module for correcting drift (within-batch correction) in metabolomics data with pooled QC samples.

### check_missingness(df: DataFrame, rows_to_check: [list](https://docs.python.org/3/library/stdtypes.html#list))

This function checks for NAs in the data frame inside some user-provided rows.

* **Parameters:**
  * **df** (*samples as rows* *,* *features as columns*)
  * **rows_to_check** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* *rows that should be checked for missingness*)
* **Returns:**
  **Boolean**
* **Return type:**
  True if there is missingness, False if there is no missingness.

### cpca_centroid(df: DataFrame, sample_rows, qc_rows: [list](https://docs.python.org/3/library/stdtypes.html#list), log_transform: [bool](https://docs.python.org/3/library/functions.html#bool) = True)

### run_cpca_drift_correction(df: DataFrame, sample_rows, qc_rows, n_comps: [int](https://docs.python.org/3/library/functions.html#int) = 1) → DataFrame

Corrects technical drift using Common Principal Components Analysis (CPCA).
Adapted from [https://github.com/m-baralt/metabolomics_incident_diabetes](https://github.com/m-baralt/metabolomics_incident_diabetes).

* **Parameters:**
  * **df** (*samples as rows* *,* *features as columns*)
  * **sample_rows** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* *sample row indices*)
  * **qc_rows** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* *QC row indices*)
  * **n_comps** (*number* *of* *common principal components to remove* *(**default 1* *)*)
* **Returns:**
  * *Full input DataFrame with corrected values applied to the intensity rows*
  *  *(sample_rows + qc_rows). Rows outside those arguments are returned unchanged.*

### run_loess_drift_correction(data, qc_rows, sample_rows, sample_order: DataFrame, filter_percent: [float](https://docs.python.org/3/library/functions.html#float) = None, qc_min_threshold: [int](https://docs.python.org/3/library/functions.html#int) = 4, always_use_default=False, default=0.75)

Perform QC-based drift correction across multiple features using
LOESS regression and spline interpolation.

For each feature:
1. Extract QC intensities and corresponding injection order.
2. Optionally filter features based on QC completeness.
3. Compute QC relative standard deviation (RSD).
4. If sufficient QC points exist, estimate a drift curve using

> qc_rlsc_loess, finding the best alpha smoothing span
> with leave-one-out cross validation (LOOCV).
1. Normalize all intensities by dividing by the drift curve and
   scaling to the QC median.
2. Record drift parameters and correction metadata.

* **Parameters:**
  * **data** ([*pandas.DataFrame*](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html#pandas.DataFrame)) – Input intensity matrix with samples as rows and features as columns.
  * **qc_rows** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* [*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Row indices corresponding to QC injections.
  * **sample_rows** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* [*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Row indices corresponding to biological samples.
  * **sample_order** ([*pandas.DataFrame*](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html#pandas.DataFrame)) – Table mapping file names to injection order. Must contain
    columns “File Name” and “Sample ID”. Sample ID must be
    numeric.
  * **filter_percent** ([*float*](https://docs.python.org/3/library/functions.html#float) *,* *optional*) – Minimum proportion of QC values that must be non-missing for
    a feature to be retained (e.g., 0.6 means at least 60% of QCs
    must be present).
  * **qc_min_threshold** ([*int*](https://docs.python.org/3/library/functions.html#int) *,* *optional*) – Minimum number of QC values required to perform drift
    correction. Features with fewer QCs are returned uncorrected.
  * **always_use_default** ([*bool*](https://docs.python.org/3/library/functions.html#bool) *,* *optional*) – If True, the alpha value 0.75 is used for the smoothing span.
    LOOCV is skipped. This option is less computationally heavy.
* **Returns:**
  * **corrected_df** (*pandas.DataFrame*) – Full input DataFrame with corrected values applied to
    sample_rows + qc_rows. Rows outside those arguments are
    returned unchanged. Features that fail QC requirements are
    returned unchanged.
  * **correction_info** (*dict*) – Dictionary keyed by feature name, containing:
    - ‘alpha’: selected LOESS alpha (or None if skipped)
    - ‘drift_curve’: the evaluated drift correction vector
    - ‘y_qc’: QC intensities used
    - ‘x_qc’: QC injection orders
    - ‘rsd_qc’: QC relative standard deviation
    - ‘median’: QC median intensity (used for scaling)
    - ‘y_all’: original intensities
    - ‘status’: “corrected”, “skipped_due_to_few_qcs”, or error note

### Notes

- Features with insufficient QC points are not corrected.
- High QC RSD (>20%) is flagged but does not prevent correction.
- Drift correction rescales intensities so that QC medians remain
  unchanged.
- Sample names in data and sample_order must match exactly.

### qc_rlsc_loess(x_qc, y_qc, x_all, always_use_default=False, default=0.75, alpha_candidates=array([0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.]))

Estimate a QC-based drift curve using LOESS smoothing with
leave-one-out cross-validation (LOOCV) to select the optimal
smoothing span (alpha).

This function:
1. Tests multiple LOESS spans (alpha values).
2. Fits LOESS to QC points for each candidate span.
3. Performs LOOCV to compute prediction error for each span.
4. Selects the alpha producing the lowest LOOCV error.
5. Fits LOESS once more using the best alpha.
6. Interpolates the LOESS fit to all sample injection orders

> using a cubic spline, with clamping outside the QC range.
* **Parameters:**
  * **x_qc** (*array-like*) – Injection order of QC samples (numeric).
  * **y_qc** (*array-like*) – Intensity values of QC samples corresponding to x_qc.
  * **x_all** (*array-like*) – Injection order of all samples (QCs + regular samples) in
    the same order as data rows.
  * **always_use_default** ([*bool*](https://docs.python.org/3/library/functions.html#bool) *,* *optional*) – If True, the alpha value 0.75 is used for all values.
    LOOCV is skipped. This option is less computationally heavy.
  * **alpha_candidates** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* [*float*](https://docs.python.org/3/library/functions.html#float) *,* *optional*) – List of LOESS smoothing parameters (fractions of data used
    in local regression) to evaluate during optimization.
    Default is [0.4, 0.6, 0.8, 1.0].
* **Returns:**
  * **drift_curve** (*ndarray*) – The estimated drift correction curve evaluated at each
    injection order in x_all. Values outside the QC range
    are clamped to the nearest in-range LOESS value.
  * **best_alpha** (*float*) – The alpha value producing the lowest LOOCV error.

### Notes

- LOESS fits enforce a minimum fraction of (λ + 1) / n, with λ=1
  for linear LOESS.
- CubicSpline is used for interpolation without extrapolation.
  Out-of-range values are manually clamped.
- Drift curve values are clipped to be strictly positive
  (minimum 1e-6) to prevent division instability.

## Submodules

## acore.drift_correction.cpca_drift_correction module

Functions for metabolomics drift correction by
Common Principal Components Analysis (CPCA).

### check_missingness(df: DataFrame, rows_to_check: [list](https://docs.python.org/3/library/stdtypes.html#list))

This function checks for NAs in the data frame inside some user-provided rows.

* **Parameters:**
  * **df** (*samples as rows* *,* *features as columns*)
  * **rows_to_check** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* *rows that should be checked for missingness*)
* **Returns:**
  **Boolean**
* **Return type:**
  True if there is missingness, False if there is no missingness.

### run_cpca_drift_correction(df: DataFrame, sample_rows, qc_rows, n_comps: [int](https://docs.python.org/3/library/functions.html#int) = 1) → DataFrame

Corrects technical drift using Common Principal Components Analysis (CPCA).
Adapted from [https://github.com/m-baralt/metabolomics_incident_diabetes](https://github.com/m-baralt/metabolomics_incident_diabetes).

* **Parameters:**
  * **df** (*samples as rows* *,* *features as columns*)
  * **sample_rows** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* *sample row indices*)
  * **qc_rows** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* *QC row indices*)
  * **n_comps** (*number* *of* *common principal components to remove* *(**default 1* *)*)
* **Returns:**
  * *Full input DataFrame with corrected values applied to the intensity rows*
  *  *(sample_rows + qc_rows). Rows outside those arguments are returned unchanged.*

### cpca_centroid(df: DataFrame, sample_rows, qc_rows: [list](https://docs.python.org/3/library/stdtypes.html#list), log_transform: [bool](https://docs.python.org/3/library/functions.html#bool) = True)

## acore.drift_correction.loess_drift_correction module

Functions for metabolomics drift correction.

### filter_features_by_qc(df: DataFrame, qc_rows: [list](https://docs.python.org/3/library/stdtypes.html#list), threshold: [float](https://docs.python.org/3/library/functions.html#float) = 0.5) → DataFrame

Filter features in a DataFrame based on quality control (QC) completeness.

This function removes columns (features) that do not meet a minimum number of valid
(non-missing) QC values. The minimum number of required valid values is computed as
ceil(n_qc \* (1 - threshold)), where n_qc is the number of QC rows.

* **Parameters:**
  * **df** ([*pandas.DataFrame*](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html#pandas.DataFrame)) – Input DataFrame with samples as rows and features as columns.
  * **qc_rows** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – List of row indices corresponding to QC samples.
  * **threshold** ([*float*](https://docs.python.org/3/library/functions.html#float)) – Fraction (between 0 and 1) indicating the maximum allowed proportion
    of missing QC values per feature. For example, threshold=0.6 allows
    up to 60% of QC values to be missing; so a feature with 3 out of 5
    QC values present (40% present, 60% missing) would be retained.
    Defaults to 0.5.
* **Returns pandas.DataFrame:**
  Filtered DataFrame containing only columns with sufficient valid QC values.
* **Raises:**
  [**ValueError**](https://docs.python.org/3/library/exceptions.html#ValueError) – If threshold is not between 0 and 1.

**Example**

> import pandas as pd
> import numpy as np

> df = pd.DataFrame({
> : ‘A’: [1.0, np.nan, 2.0],
>   ‘B’: [np.nan, 3.0, np.nan],
>   ‘C’: [4.0, np.nan, 6.0]

> }, index=[‘QC1’, ‘QC2’, ‘QC3’])

> filtered = filter_features_by_qc(df, qc_rows=[‘QC1’, ‘QC2’, ‘QC3’], threshold=0.5)
> print(filtered)
> # Output: columns with at least 2 valid QC values (since ceil(3 \* (1 - 0.5)) = 2)

### qc_rlsc_loess(x_qc, y_qc, x_all, always_use_default=False, default=0.75, alpha_candidates=array([0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.]))

Estimate a QC-based drift curve using LOESS smoothing with
leave-one-out cross-validation (LOOCV) to select the optimal
smoothing span (alpha).

This function:
1. Tests multiple LOESS spans (alpha values).
2. Fits LOESS to QC points for each candidate span.
3. Performs LOOCV to compute prediction error for each span.
4. Selects the alpha producing the lowest LOOCV error.
5. Fits LOESS once more using the best alpha.
6. Interpolates the LOESS fit to all sample injection orders

> using a cubic spline, with clamping outside the QC range.
* **Parameters:**
  * **x_qc** (*array-like*) – Injection order of QC samples (numeric).
  * **y_qc** (*array-like*) – Intensity values of QC samples corresponding to x_qc.
  * **x_all** (*array-like*) – Injection order of all samples (QCs + regular samples) in
    the same order as data rows.
  * **always_use_default** ([*bool*](https://docs.python.org/3/library/functions.html#bool) *,* *optional*) – If True, the alpha value 0.75 is used for all values.
    LOOCV is skipped. This option is less computationally heavy.
  * **alpha_candidates** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* [*float*](https://docs.python.org/3/library/functions.html#float) *,* *optional*) – List of LOESS smoothing parameters (fractions of data used
    in local regression) to evaluate during optimization.
    Default is [0.4, 0.6, 0.8, 1.0].
* **Returns:**
  * **drift_curve** (*ndarray*) – The estimated drift correction curve evaluated at each
    injection order in x_all. Values outside the QC range
    are clamped to the nearest in-range LOESS value.
  * **best_alpha** (*float*) – The alpha value producing the lowest LOOCV error.

### Notes

- LOESS fits enforce a minimum fraction of (λ + 1) / n, with λ=1
  for linear LOESS.
- CubicSpline is used for interpolation without extrapolation.
  Out-of-range values are manually clamped.
- Drift curve values are clipped to be strictly positive
  (minimum 1e-6) to prevent division instability.

### run_loess_drift_correction(data, qc_rows, sample_rows, sample_order: DataFrame, filter_percent: [float](https://docs.python.org/3/library/functions.html#float) = None, qc_min_threshold: [int](https://docs.python.org/3/library/functions.html#int) = 4, always_use_default=False, default=0.75)

Perform QC-based drift correction across multiple features using
LOESS regression and spline interpolation.

For each feature:
1. Extract QC intensities and corresponding injection order.
2. Optionally filter features based on QC completeness.
3. Compute QC relative standard deviation (RSD).
4. If sufficient QC points exist, estimate a drift curve using

> qc_rlsc_loess, finding the best alpha smoothing span
> with leave-one-out cross validation (LOOCV).
1. Normalize all intensities by dividing by the drift curve and
   scaling to the QC median.
2. Record drift parameters and correction metadata.

* **Parameters:**
  * **data** ([*pandas.DataFrame*](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html#pandas.DataFrame)) – Input intensity matrix with samples as rows and features as columns.
  * **qc_rows** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* [*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Row indices corresponding to QC injections.
  * **sample_rows** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* [*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Row indices corresponding to biological samples.
  * **sample_order** ([*pandas.DataFrame*](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html#pandas.DataFrame)) – Table mapping file names to injection order. Must contain
    columns “File Name” and “Sample ID”. Sample ID must be
    numeric.
  * **filter_percent** ([*float*](https://docs.python.org/3/library/functions.html#float) *,* *optional*) – Minimum proportion of QC values that must be non-missing for
    a feature to be retained (e.g., 0.6 means at least 60% of QCs
    must be present).
  * **qc_min_threshold** ([*int*](https://docs.python.org/3/library/functions.html#int) *,* *optional*) – Minimum number of QC values required to perform drift
    correction. Features with fewer QCs are returned uncorrected.
  * **always_use_default** ([*bool*](https://docs.python.org/3/library/functions.html#bool) *,* *optional*) – If True, the alpha value 0.75 is used for the smoothing span.
    LOOCV is skipped. This option is less computationally heavy.
* **Returns:**
  * **corrected_df** (*pandas.DataFrame*) – Full input DataFrame with corrected values applied to
    sample_rows + qc_rows. Rows outside those arguments are
    returned unchanged. Features that fail QC requirements are
    returned unchanged.
  * **correction_info** (*dict*) – Dictionary keyed by feature name, containing:
    - ‘alpha’: selected LOESS alpha (or None if skipped)
    - ‘drift_curve’: the evaluated drift correction vector
    - ‘y_qc’: QC intensities used
    - ‘x_qc’: QC injection orders
    - ‘rsd_qc’: QC relative standard deviation
    - ‘median’: QC median intensity (used for scaling)
    - ‘y_all’: original intensities
    - ‘status’: “corrected”, “skipped_due_to_few_qcs”, or error note

### Notes

- Features with insufficient QC points are not corrected.
- High QC RSD (>20%) is flagged but does not prevent correction.
- Drift correction rescales intensities so that QC medians remain
  unchanged.
- Sample names in data and sample_order must match exactly.
