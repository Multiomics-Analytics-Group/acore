"""
Functions for metabolomics drift correction.
"""

import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import CubicSpline
import numpy as np


def filter_features_by_qc(
    df: pd.DataFrame, qc_cols: list, threshold: float = 0.5
) -> pd.DataFrame:
    """
    Filter features in a DataFrame based on quality control (QC) completeness.

    This function removes rows (features) that do not meet a minimum number of valid
    (non-missing) QC values. The minimum number of required valid values is computed as
    `ceil(n_qc * (1 - threshold))`, where `n_qc` is the number of QC columns.

    :param pandas.DataFrame df:
        Input DataFrame containing feature data and QC columns.
    :param list qc_cols:
        List of column names corresponding to QC measurements.
    :param float threshold:
        Fraction (between 0 and 1) indicating the maximum allowed proportion
        of missing QC values per feature. For example, `threshold=0.5` allows
        up to 50% missing QC values. Defaults to 0.5.

    :returns pandas.DataFrame:
        Filtered DataFrame containing only rows with sufficient valid QC values.

    :raises ValueError:
        If `threshold` is not between 0 and 1.

    **Example**

        import pandas as pd
        import numpy as np

        df = pd.DataFrame({
            'feature': ['A', 'B', 'C'],
            'QC1': [1.0, np.nan, 2.0],
            'QC2': [np.nan, 3.0, np.nan],
            'QC3': [4.0, np.nan, 6.0]
        })

        filtered = filter_features_by_qc(df, qc_cols=['QC1', 'QC2', 'QC3'], threshold=0.5)
        print(filtered)
        # Output: rows with at least 2 valid QC values (since ceil(3 * (1 - 0.5)) = 2)

    """
    n_qc = len(qc_cols)
    min_valid = int(np.ceil(n_qc * (1 - threshold)))  # minimum valid QC points
    valid_counts = df[qc_cols].notna().sum(axis=1)

    return df.loc[valid_counts >= min_valid]


def qc_rlsc_loess(
    x_qc,
    y_qc,
    x_all,
    use_default=False,
    alpha_candidates=np.arange(0.4, 1.01, 0.05),
    print_logs=False,
):
    """
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
       using a cubic spline, with clamping outside the QC range.

    Parameters
    ----------
    x_qc : array-like
        Injection order of QC samples (numeric).
    y_qc : array-like
        Intensity values of QC samples corresponding to `x_qc`.
    x_all : array-like
        Injection order of all samples (QCs + regular samples) in
        the same order as data columns.
    use_default: bool, optional
        If True, the alpha value 0.75 is used for all values.
        LOOCV is skipped. This option is less computationally heavy.
    alpha_candidates : list of float, optional
        List of LOESS smoothing parameters (fractions of data used
        in local regression) to evaluate during optimization.
        Default is [0.4, 0.6, 0.8, 1.0].
    print_logs : bool, optional
        If True, prints diagnostic information during optimization.

    Returns
    -------
    drift_curve : ndarray
        The estimated drift correction curve evaluated at each
        injection order in `x_all`. Values outside the QC range
        are clamped to the nearest in-range LOESS value.
    best_alpha : float
        The alpha value producing the lowest LOOCV error.

    Notes
    -----
    - LOESS fits enforce a minimum fraction of (λ + 1) / n, with λ=1
      for linear LOESS.
    - CubicSpline is used for interpolation without extrapolation.
      Out-of-range values are manually clamped.
    - Drift curve values are clipped to be strictly positive
      (minimum 1e-6) to prevent division instability.
    """

    best_alpha = None
    best_loocv_error = np.inf
    n = len(x_qc)

    if not use_default:  # do leave one out cross validation for determining best alpha
        for alpha in alpha_candidates:
            span = max(
                (1 + 1) / n, alpha
            )  # minimum alpha is (λ+1)/n, λ=1 for linear LOESS
            # Compute LOESS on full QC data
            loess_fit = lowess(y_qc, x_qc, frac=span, it=0, return_sorted=False)

            # Leave-one-out cross validation
            errors = []
            for i in range(n):
                x_cv = np.delete(x_qc, i)
                y_cv = np.delete(y_qc, i)
                loess_cv = lowess(
                    y_cv, x_cv, frac=span, it=0, return_sorted=False, xvals=[x_qc[i]]
                )
                errors.append((y_qc[i] - loess_cv[0]) ** 2)
            loocv_error = np.mean(errors)

            if loocv_error < best_loocv_error:
                best_loocv_error = loocv_error
                best_alpha = alpha

    # Fit LOESS again with best alpha for final curve
    if best_alpha is None or use_default:
        # fallback to reasonable default if optimization failed
        best_alpha = 0.75
        if print_logs and not use_default:
            print(
                f"Warning: LOESS optimization failed for n={n}. Using default alpha=0.75."
            )

    span = max((1 + 1) / n, best_alpha)
    loess_fit = lowess(y_qc, x_qc, frac=span, it=0, return_sorted=False)

    # Cubic spline interpolation for all points (samples + QCs)
    # Restrict to interpolation range only
    cs = CubicSpline(
        x_qc, loess_fit, extrapolate=False
    )  # No extrapolation outside QC range
    drift_curve = cs(x_all)

    # Optional: Clip drift_curve to the edge values to prevent NaNs or negatives
    x_min, x_max = np.min(x_qc), np.max(x_qc)

    # For values outside QC range, hold the first/last fitted value (clamping)
    drift_curve[x_all < x_min] = loess_fit[0]
    drift_curve[x_all > x_max] = loess_fit[-1]
    # Final safeguard: ensure no negatives
    drift_curve = np.clip(drift_curve, a_min=1e-6, a_max=None)

    return drift_curve, best_alpha


def run_drift_correction(
    data,
    qc_cols,
    sample_cols,
    sample_order: pd.DataFrame,
    feature_name_col: str = None,
    filter_percent: float = None,
    qc_min_threshold: int = 4,
    print_logs=False,
    use_default=False,
):
    """
    Perform QC-based drift correction across multiple features using
    LOESS regression and spline interpolation.

    For each feature:
    1. Extract QC intensities and corresponding injection order.
    2. Optionally filter features based on QC completeness.
    3. Compute QC relative standard deviation (RSD).
    4. If sufficient QC points exist, estimate a drift curve using
       `qc_rlsc_loess`, finding the best alpha smoothing span
       with leave-one-out cross validation (LOOCV).
    5. Normalize all intensities by dividing by the drift curve and
       scaling to the QC median.
    6. Record drift parameters and correction metadata.

    Parameters
    ----------
    data : pandas.DataFrame
        Input intensity matrix, with features as rows and sample/QC
        names as columns.
    qc_cols : list of str
        Names of columns corresponding to QC injections.
    sample_cols : list of str
        Names of columns corresponding to biological samples.
    sample_order : pandas.DataFrame
        Table mapping file names to injection order. Must contain
        columns "File Name" and "Sample ID". Sample ID must be
        numeric.
    feature_name_col : str, optional
        Name of the column in `data` containing feature identifiers.
        If None, the index is used.
    filter_percent : float, optional
        Minimum proportion of QC values that must be non-missing for
        a feature to be retained (e.g., 0.5 means at least 50% of QCs
        must be present).
    qc_min_threshold : int, optional
        Minimum number of QC values required to perform drift
        correction. Features with fewer QCs are returned uncorrected.
    print_logs : bool, optional
        If True, prints status messages, RSD warnings, and errors.
    use_default: bool, optional
        If True, the alpha value 0.75 is used for the smoothing span.
        LOOCV is skipped. This option is less computationally heavy.

    Returns
    -------
    corrected_df : pandas.DataFrame
        A feature x sample matrix of drift-corrected intensities.
        All sample and QC columns are included. Features that fail
        QC requirements are returned unchanged.
    correction_info : dict
        Dictionary keyed by feature name, containing:
        - 'alpha': selected LOESS alpha (or None if skipped)
        - 'drift_curve': the evaluated drift correction vector
        - 'y_qc': QC intensities used
        - 'x_qc': QC injection orders
        - 'rsd_qc': QC relative standard deviation
        - 'median': QC median intensity (used for scaling)
        - 'y_all': original intensities
        - 'status': "corrected", "skipped_due_to_few_qcs", or error note

    Notes
    -----
    - Features with insufficient QC points are not corrected.
    - High QC RSD (>20%) is flagged but does not prevent correction.
    - Drift correction rescales intensities so that QC medians remain
      unchanged.
    - Sample names in `data` and `sample_order` must match exactly.
    """
    df = data.copy()
    correction_info = {}  # Feature name -> dict with 'alpha' and 'drift_curve'
    all_cols = sample_cols + qc_cols
    corrected_df = pd.DataFrame(
        index=df.index, columns=all_cols, dtype=float
    )  # Initialises new df to be filled (for now has all nan), all columns including QC should be corrected
    if feature_name_col:
        df["TempName"] = df[
            feature_name_col
        ]  # naming it TempName so it does not interfere with any other col called Name
    else:
        df["TempName"] = df.index

    try:
        sample_order["Sample ID"] = pd.to_numeric(
            sample_order["Sample ID"], errors="raise"
        ).astype(int)
    except ValueError as e:
        raise ValueError(
            "Sample ID column in sample order data contains non-integer values. All values in this column need to be numbers."
        ) from e

    # Injection order mapping (sample name -> injection order)
    injection_order_map = dict(
        zip(sample_order["File Name"], sample_order["Sample ID"])
    )

    x_all = np.array(
        [injection_order_map.get(sample, np.nan) for sample in all_cols]
    )  # Order of samples in order of file names
    print(x_all)
    if np.isnan(x_all).any():
        print(
            "Warning: some of your samples don't have an associated sample order. They will be skipped.",
            "Make sure the names of your samples are identical in both your data frame and your sample order data.",
        )

    if filter_percent is not None:  # Filter features by QC completeness
        df = filter_features_by_qc(
            df, qc_cols, threshold=filter_percent
        )  # Only keeping features that have min half of all QCs not nan

    # Loop through all features in filtered df
    for feature_idx, row in df[all_cols].iterrows():

        # Y VALUES ARE INTENSITIES
        # X VALUES ARE RUN ORDER IDX

        cid_value = df.at[feature_idx, "TempName"]  # feature ID

        # Intensities (y) for all cols in the samples (t, c and nan)
        y_all = row.values.astype(float)

        # Get intensities (y) and run order idx (x) for all QC data points
        y_qc = row[qc_cols].values.astype(
            float
        )  # Gets all intensity values for row for all qcs
        x_qc = np.array(
            [injection_order_map.get(s, np.nan) for s in qc_cols]
        )  # Gets the order for those

        # Remove NaNs from QC for fitting
        valid_mask = ~np.isnan(y_qc) & ~np.isnan(
            x_qc
        )  # Figure out which ones are nan (mask, has True/False in it)
        y_qc_valid = y_qc[valid_mask]  # Same as y_qc but only non-nans
        x_qc_valid = x_qc[valid_mask]

        # Calculate RSD
        rsd_qc = 100 * np.std(y_qc_valid) / np.mean(y_qc_valid)
        if rsd_qc > 30 and print_logs:
            print(
                f"Flagging feature {cid_value} due to too high QC RSD. (Feature not skipped)"
            )  # RSD>20% suggests qc intensities reflect instrument drift, not biology
            print("RSD:", rsd_qc)

        if len(y_qc_valid) < qc_min_threshold:
            # Skip correction but preserve the row
            if print_logs:
                print("skipping correction due to few QCs:", cid_value)
            corrected_df.loc[feature_idx] = y_all  # Insert the uncorrected values
            correction_info[cid_value] = {
                "alpha": None,
                "drift_curve": None,
                "y_qc": y_qc_valid.tolist(),
                "x_qc": x_qc_valid.tolist(),
                "rsd_qc": rsd_qc,
                "status": "skipped_due_to_few_qcs",
            }
            continue

        try:  # DRIFT CORRECTION
            # Calculate drift curve
            drift_curve, best_alpha = qc_rlsc_loess(
                x_qc_valid,
                y_qc_valid,
                x_all,
                print_logs=print_logs,
                use_default=use_default,
            )  # Calculate curve and alpha value
            median_qc = np.median(y_qc_valid)
            # print("Any negatives in drift_curve?", np.any(drift_curve < 0))
            # print("median", median_qc)

            # Remove NaNs from y_all for indexing (just keep for drift correction) -> only if they are not nan in either data or drift curve
            valid_mask = ~np.isnan(y_all) & ~np.isnan(drift_curve)

            # Normalize all intensities using drift curve
            corrected = np.full_like(
                y_all, np.nan, dtype=float
            )  # Initialise corrected array
            corrected[valid_mask] = (
                y_all[valid_mask] / drift_curve[valid_mask]
            ) * median_qc  # The ones that were nan don't get corrected

            corrected_df.loc[feature_idx] = corrected

            correction_info[cid_value] = {
                "alpha": best_alpha,
                "drift_curve": drift_curve.tolist(),  # Convert to list for JSON/pickle
                "y_qc": y_qc_valid.tolist(),  # Store QC values
                "x_qc": x_qc_valid.tolist(),  # Store QC injection orders
                "rsd_qc": rsd_qc,
                "median": median_qc,
                "y_all": y_all,
                "status": "corrected",
            }

            if print_logs:
                print(f"Corrected {cid_value} with alpha {best_alpha}.")

        except Exception as e:
            print(f"Skipping feature {cid_value} due to error: {e}")

            continue

    corrected_df["Name"] = df["TempName"]
    return corrected_df, correction_info
