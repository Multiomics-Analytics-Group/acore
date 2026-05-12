# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: acore-dev
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Metabolomics drift correction
#

# %% [markdown]
# Within-batch correction of metabolomics data.
# This script shows how to correct for instrumental drift based on pooled QC data samples.

# %% [markdown]
# ## 1. LOESS smoothing-based drift correction

# %%
# #%pip install acore

# %%
import importlib
import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import acore
from acore import drift_correction as dc

importlib.reload(acore)


# %%
def plot_loess_example_curve(
    df: pd.DataFrame,
    feature_idx: int,
    sample_cols: list,
    qc_cols: list,
    sample_order: pd.DataFrame,
    feature_name_col: str = None,
    show_corrected: bool = True,
    alpha: float = None,  # fixed smoothing span; if None, selected by LOOCV
):
    """
    Plot the raw intensities, LOESS drift curve, and optionally the corrected
    intensities for a single feature according to the loess drift correction function.
    Useful for inspecting drift behaviour before running full drift correction.

    The drift curve is estimated with the same method used by
    ldc.run_drift_correction: LOESS is fit to the QC points and then
    interpolated across all injection positions via a cubic spline
    (qc_rlsc_loess).

    Parameters
    ----------
    df : pd.DataFrame
        Feature matrix with features as rows and samples/QCs as columns.
    feature_idx : int
        Row index of the feature to plot.
    sample_cols : list of str
        Column names of the biological samples.
    qc_cols : list of str
        Column names of the pooled QC samples.
    sample_order : pd.DataFrame
        Injection-order table with columns "File Name" and "Sample ID"
        (integer run order).
    feature_name_col : str, optional
        Column in df containing feature identifiers used in the plot title.
        If None, the row index is used.
    show_corrected : bool, optional
        If True (default), overlays drift-corrected sample intensities as
        diamond markers.
    alpha : float, optional
        LOESS smoothing span (0 < α ≤ 1). If None (default), the optimal span
        is selected automatically by leave-one-out cross-validation over
        α ∈ [0.40, 1.00]. The selected value is shown in the legend.
    """

    feature_row = df.iloc[feature_idx]
    all_cols = sample_cols + qc_cols

    order_dict = sample_order.set_index("File Name")["Sample ID"].to_dict()
    x_all = np.array([order_dict.get(c, np.nan) for c in all_cols])
    y_all = feature_row[all_cols].astype(float).values

    n_s = len(sample_cols)
    x_sample, y_sample = x_all[:n_s], y_all[:n_s]
    x_qc_arr, y_qc_arr = x_all[n_s:], y_all[n_s:]

    valid_sample = ~np.isnan(x_sample) & ~np.isnan(y_sample)
    valid_qc = ~np.isnan(x_qc_arr) & ~np.isnan(y_qc_arr)
    x_s_v, y_s_v = x_sample[valid_sample], y_sample[valid_sample]
    x_qc_v, y_qc_v = x_qc_arr[valid_qc], y_qc_arr[valid_qc]

    if feature_name_col and feature_name_col in df.columns:
        feature_name = feature_row[feature_name_col]
    else:
        feature_name = f"index {feature_idx}"

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.scatter(x_s_v, y_s_v, label="Samples", color="steelblue", alpha=0.7, zorder=3)
    ax.scatter(
        x_qc_v, y_qc_v, label="QC", color="firebrick", edgecolor="k", s=60, zorder=4
    )

    if len(x_qc_v) >= 4:
        if alpha is not None:
            # Pass a single-element candidate list to skip LOOCV and use the given alpha directly
            drift_curve, best_alpha = dc.qc_rlsc_loess(
                x_qc_v, y_qc_v, x_all, use_default=True, default=alpha
            )
        else:
            drift_curve, best_alpha = dc.qc_rlsc_loess(x_qc_v, y_qc_v, x_all)

        valid_curve = ~np.isnan(x_all) & ~np.isnan(drift_curve)
        sort_idx = np.argsort(x_all[valid_curve])
        ax.plot(
            x_all[valid_curve][sort_idx],
            drift_curve[valid_curve][sort_idx],
            label=f"LOESS drift curve (α={best_alpha:.2f})",
            color="black",
            lw=2,
            zorder=5,
        )

        if show_corrected:
            median_qc = np.median(y_qc_v)
            drift_at_samples = drift_curve[:n_s][valid_sample]
            corrected = (y_s_v / drift_at_samples) * median_qc
            ax.scatter(
                x_s_v,
                corrected,
                label="Corrected samples",
                color="lightsteelblue",
                marker="D",
                s=40,
                alpha=0.9,
                zorder=3,
            )
    else:
        print(f"Not enough valid QC points for LOESS ({len(x_qc_v)} found, need ≥4).")

    ax.set_xlabel("Injection Order")
    ax.set_ylabel("Intensity")
    ax.set_title(f"Drift correction example ({feature_name})")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


# %%
def pca_for_cpca_drift(
    df: pd.DataFrame,
    sample_cols,  # list of col names, OR dict {group_name: [col names]}
    qc_cols: list,
    log_transform: bool = True,
    title: str = "PCA",
):
    """
    PCA of samples and QC samples.

    Parameters
    ----------
    df          : feature matrix with features as rows, samples as columns
    sample_cols : list of sample column names (all one group),
                  or dict {group_name: [col names]} for multiple groups
    qc_cols     : list of QC column names
    log_transform : apply log1p before scaling
    title       : plot title
    """
    # Normalise sample_cols to a dict
    if isinstance(sample_cols, list):
        sample_groups = {"Samples": sample_cols}
    else:
        sample_groups = sample_cols

    # Build ordered column + label arrays
    all_cols, labels = [], []
    for group, cols in sample_groups.items():
        for c in cols:
            if c in df.columns:
                all_cols.append(c)
                labels.append(group)
    for c in qc_cols:
        if c in df.columns:
            all_cols.append(c)
            labels.append("QC")

    # Coerce to numeric, drop features with any NaN
    X = df[all_cols].apply(pd.to_numeric, errors="coerce").dropna(axis=0)

    if log_transform:
        X = np.log1p(X.clip(lower=0))

    X_scaled = StandardScaler().fit_transform(X.T.values.astype(float))

    pca = PCA(n_components=2)
    coords = pca.fit_transform(X_scaled)
    pc1_var = pca.explained_variance_ratio_[0] * 100
    pc2_var = pca.explained_variance_ratio_[1] * 100

    palette = plt.cm.tab10.colors
    color_map = {g: palette[i % len(palette)] for i, g in enumerate(sample_groups)}
    color_map["QC"] = "black"

    fig, ax = plt.subplots(figsize=(8, 6))
    for group in list(sample_groups) + ["QC"]:
        mask = [_label == group for _label in labels]
        pts = coords[mask]
        is_qc = group == "QC"
        ax.scatter(
            pts[:, 0],
            pts[:, 1],
            label=group,
            color=color_map[group],
            marker="D" if is_qc else "o",
            s=70 if is_qc else 50,
            edgecolors="k" if is_qc else "none",
            alpha=0.9 if is_qc else 0.75,
            zorder=5 if is_qc else 3,
        )

    ax.set_xlabel(f"PC1 ({pc1_var:.1f}%)")
    ax.set_ylabel(f"PC2 ({pc2_var:.1f}%)")
    ax.set_title(title)
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.axvline(0, color="grey", lw=0.5, ls="--")
    ax.legend(framealpha=0.8)
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.show()


# %% [markdown]
# ### Load in data

# %%
df = pd.read_excel(
    "../../example_data/aradopsis_seedling_lipids/kehelpannala_AnArabidopsisLipid_seedling.xlsx"
)
sample_order = pd.read_csv(
    "../../example_data/aradopsis_seedling_lipids/seedling_art_sample_order.csv"
)

# %%
df

# %% [markdown]
# We also have (artificial) data that has the order in which our samples were run. This information is crucial for the drift correction algorithm.

# %%
sample_order.sort_values("Sample ID")

# %% [markdown]
# ### Run drift correction with LOESS smoothing

# %% [markdown]
# We can now correct our data for experimental drift.
#
# With the acore LOESS drift-correction function, a LOESS (locally estimated regression) smoother is applied separately to the features in the data to model slow temporal trends, and the resulting smooth trend is used to correct the data.
#
# Before the function estimation and correction, the data can filtered, to remove features that have too many missing values in the QC samples.

# %%
# First, we can create a dictionary for our sample names, ordering them into groups, to make the upcoming function call easier.

data_groups = {
    "SD1": ["Sd1_1", "Sd1_2", "Sd1_3", "Sd1_4", "Sd1_5", "Sd1_6"],
    "SD2": [
        "Sd2_1",
        "Sd2_2",
        "Sd2_3",
        "Sd2_4",
        "Sd2_5",
        "Sd2_6",
        "Sd2_7",
        "Sd2_8",
        "Sd2_9",
    ],
    "SD3": [
        "Sd3-1",
        "Sd3-2",
        "Sd3-3",
        "Sd3-4",
        "Sd3-5",
        "Sd3-6",
        "Sd3-7",
        "Sd3-8",
        "Sd3-9",
    ],
    "QC": [
        "PBQC_Sd_1",
        "PBQC_Sd_2",
        "PBQC_Sd_3",
        "PBQC_Sd_4",
        "PBQC_Sd_5",
        "PBQC_Sd_6",
    ],
}

# %%
# Separate groups
sample_cols = data_groups["SD1"] + data_groups["SD2"] + data_groups["SD3"]
qc_cols = data_groups["QC"]

# Create a sub-df consisting only of the interesting columns, omitting metadata.
df_dc = df[sample_cols + qc_cols]

# %% [markdown]
# Now we can run the drift correction.

# %%
corrected_df, correction_info = dc.run_loess_drift_correction(
    df,
    qc_cols,
    sample_cols,
    sample_order=sample_order,
    feature_name_col=None,
    filter_percent=0.5,
)

### Explanation of the parameters chosen
# - feature_name_col = the name of the column containing feature names, if there is one.
#   This information is used for logging and showing outputs, it's not required for the functioning of the method.
#   Here, there is no feature name column available, so "None" is used.
# - filter_percent =  the minimum percentage of values that must be present for this feature to be retained.
#   If the percentage of non-missing is below this, the feature will be filtered out.
#   If this parameter is set to "None", no filtering will be done.
# - print_logs = whether there should be an output for the logs of the function. Like verbose.

# %% [markdown]
# Now we can inspect our results. First, the corrected output dataframe.

# %%
corrected_df

# %%
df

# %% [markdown]
# We can also look further into the correction_info object, to see the parameters that were chosen for each feature.
#
# For example, let's check the parameters used for the 200th feature.

# %%
correction_info[200]

# %% [markdown]
# ### Plot an example curve for one feature

# %% [markdown]
# We can also plot an example feature, to see how the values have changed and what the LOESS curve would look like for the data of this feature.

# %%
plot_loess_example_curve(
    df=df,
    feature_idx=5,
    sample_cols=sample_cols,
    qc_cols=qc_cols,
    sample_order=sample_order,
)

# %% [markdown]
# In this plot, we can see all of the data points of this feature, ordered by the injection time. The red points are our QC samples, so we can see whether there was any instrumental drift over time. The LOESS curve is also calculated, with the smoothing value alpha chosen with leave-one-out cross validation, just like in the run_drift_correction function.
#
# Alternative values for the smoothing parameter alpha can be tested by adding the argument "alpha" and choosing a value, just like in the example below.

# %%
plot_loess_example_curve(
    df=df,
    feature_idx=5,
    sample_cols=sample_cols,
    qc_cols=qc_cols,
    sample_order=sample_order,
    alpha=0.6,
)

# %% [markdown]
# ## 2. Common Principal Components Analysis-based drift correction

# %% [markdown]
# We can also use another method for correcting drift which is based on Common Principal Components Analysis (CPCA).
# This method is based on common principal components in defined groups of the data. It assumes that when calculating common principal components of QC samples, the drift contribution can be identified as the direction capturing maximum variance that simultaneously diagonalizes the covariance matrices of a set of classes.
#
# Therefore, the variability in the identified direction can be explained as caused by experimental drift and subtracted from all samples.
#
# Let's use different example data for demonstrating this method.

# %% [markdown]
# ### Load in data

# %%
# Load data
df = pd.read_csv(
    "../../example_data/DidacMauricio_hilic/DM_FIS2018_Hilic_pos_results2023_filled_imputed.csv"
)

# Define sample columns and qc columns
collist = list(df.columns.values)
sample_cols = []
qc_cols = []
for col in collist:
    if col.startswith("AAA"):
        sample_cols.append(col)
    elif col.startswith("QC"):
        qc_cols.append(col)

# %%
df

# %% [markdown]
# Seeing as this method is based on calculating principal components, as with PCA, there must not be any missing data.
#
# We will therefore first calculate missingness (NAs).
#
# We need to check for missing values in both the sample columns and the QC columns.

# %%
if dc.check_missingness(df, sample_cols + qc_cols):
    print(
        "There are missing values. Consider imputing first, or use LOESS drift correction instead."
    )
else:
    print("\nThere is no missingness. We can proceed with the CPCA drift correction.")

# %% [markdown]
# ### Visualise non-corrected data
#
# Now that we know we can proceed, let's visualise our data before drift correction with a PCA.

# %%
pca_for_cpca_drift(
    df,
    sample_cols,  # list of col names, OR dict {group_name: [col names]}
    qc_cols,
    log_transform=True,
    title="PCA",
)

# %% [markdown]
# We see some clear indication of instrumental drift in the QC samples.

# %% [markdown]
# ### Run drift correction based on CPCA

# %%
df_corrected = dc.run_cpca_drift_correction(df, sample_cols, qc_cols, n_comps=1)

# %% [markdown]
# Let's plot the corrected data.

# %%
pca_for_cpca_drift(  # This function is defined in the beginning of this notebook.
    df_corrected,
    sample_cols,  # list of col names, OR dict {group_name: [col names]}
    qc_cols,
    log_transform=True,
    title="PCA",
)

# %% [markdown]
# There is some change, but still a significant amount of drift is clearly visible from the PCA plot.
#
# We can play around with the n_comps variable which decides the number of components.

# %%
df_corrected_2comps = dc.run_cpca_drift_correction(df, sample_cols, qc_cols, n_comps=2)

df_corrected_3comps = dc.run_cpca_drift_correction(df, sample_cols, qc_cols, n_comps=3)

df_corrected_4comps = dc.run_cpca_drift_correction(df, sample_cols, qc_cols, n_comps=4)

# %%
pca_for_cpca_drift(
    df_corrected_2comps,
    sample_cols,  # list of col names, OR dict {group_name: [col names]}
    qc_cols,
    log_transform=True,
    title="PCA with 2 components",
)

pca_for_cpca_drift(
    df_corrected_3comps,
    sample_cols,  # list of col names, OR dict {group_name: [col names]}
    qc_cols,
    log_transform=True,
    title="PCA with 3 components",
)

pca_for_cpca_drift(
    df_corrected_4comps,
    sample_cols,  # list of col names, OR dict {group_name: [col names]}
    qc_cols,
    log_transform=True,
    title="PCA with 4 components",
)

# %% [markdown]
# The correction methods with three and four components are already looking better. Let's calculate the centroids of the QC principal components and the distance of the QC points to them, to objectively decide which number of n_comps is most favourable.

# %%
dc.cpca_centroid(df_corrected, sample_cols, qc_cols, log_transform=True)
dc.cpca_centroid(df_corrected_2comps, sample_cols, qc_cols, log_transform=True)
dc.cpca_centroid(df_corrected_3comps, sample_cols, qc_cols, log_transform=True)
dc.cpca_centroid(df_corrected_4comps, sample_cols, qc_cols, log_transform=True)

# %% [markdown]
# According to this, the CPCA method with three principal components is most favourable in this case.
#
# We can go ahead and continue our metabolomics data analysis with this data set.

# %%
df_corrected_3comps
