"""
Functions for metabolomics drift correction by
Common Principal Components Analysis (CPCA).
"""

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def check_missingness(df: pd.DataFrame, cols_to_check: list):
    """
    This function checks for NAs in the data frame inside some user-provided columns.

    Parameters
    ---------
    df      : features as rows, samples as columns
    cols_to_check    : list of columns that should be checked for missingness

    Returns
    --------
    Boolean: True if there is no missingness, False if there is missingness.
    """
    na_counts = df[cols_to_check].isna().sum()
    na_features = df[cols_to_check].isna().any(axis=1).sum()

    print(f"Features (rows) with at least one NA: {na_features} / {len(df)}")
    print(f"Samples (cols) with at least one NA:")
    print(na_counts[na_counts > 0])

    if na_counts.any():
        return True
    else:
        return False


def run_cpca_drift_correction(
    df: pd.DataFrame, sample_cols, qc_cols, n_comps: int = 1
) -> pd.DataFrame:
    """
    Corrects technical drift using Common Principal Components Analysis (CPCA).
    Adapted from https://github.com/m-baralt/metabolomics_incident_diabetes.

    Parameters
    ----------
    df          : features as rows, samples as columns
    sample_cols : list of sample column names
    qc_cols     : list of QC column names
    n_comps     : number of common principal components to remove (default 1)

    Returns
    -------
    Full drift-corrected DataFrame with metadata columns preserved.
    """
    intensity_cols = sample_cols + qc_cols
    meta_cols = [c for c in df.columns if c not in intensity_cols]

    X = df[intensity_cols].values.astype(float)

    if np.isnan(X).any():
        raise ValueError("NA values present in dataset. Consider imputation first.")

    sample_idx = list(range(len(sample_cols)))
    qc_idx = list(range(len(sample_cols), len(intensity_cols)))

    scaler = StandardScaler()
    Xs = scaler.fit_transform(X.T)  # shape: (samples, features)

    cov_sample = np.cov(Xs[sample_idx], rowvar=False)
    cov_qc = np.cov(Xs[qc_idx], rowvar=False)
    avg_cov = (cov_sample + cov_qc) / 2

    k = max(n_comps, 3)
    eigenvalues, eigenvectors = np.linalg.eigh(avg_cov)
    eigenvectors = eigenvectors[:, np.argsort(eigenvalues)[::-1]]
    cpcs = eigenvectors[:, :k]

    var_projected = np.sum((Xs @ cpcs) ** 2, axis=0)
    var_cpc = np.round(var_projected / np.sum(Xs**2), 3)[:n_comps]
    print(
        "CPC explained variance:",
        dict(zip([f"CPC{i+1}" for i in range(n_comps)], var_cpc)),
    )

    W = cpcs[:, :n_comps]
    Xs_corrected = Xs - Xs @ W @ W.T
    X_corrected = scaler.inverse_transform(
        Xs_corrected
    ).T  # back to (features, samples)

    df_corrected = pd.DataFrame(X_corrected, index=df.index, columns=intensity_cols)
    return pd.concat([df[meta_cols], df_corrected], axis=1)


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
        mask = [l == group for l in labels]
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


def cpca_centroid(
    df: pd.DataFrame,
    sample_cols,  # list of col names, OR dict {group_name: [col names]}
    qc_cols: list,
    log_transform: bool = True,
):
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

    # define qc points in coords
    for group in list(sample_groups) + ["QC"]:
        mask = [l == group for l in labels]
        qc_points = coords[mask]

    # Calculate centroid
    length = qc_points.shape[0]
    sum_x = np.sum(qc_points[:, 0])
    sum_y = np.sum(qc_points[:, 1])
    centroid = sum_x / length, sum_y / length

    distances = [math.dist(point, centroid) for point in qc_points]
    total_distance = sum(distances)
    avg_distance = total_distance / length
    print(f"The average distance of the points to the centroid is {avg_distance:.3f}.")

    return avg_distance
