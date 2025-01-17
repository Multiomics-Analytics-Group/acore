"""Decompositon plots like pca, umap, tsne, etc."""

from typing import Optional

import matplotlib
import pandas as pd
import sklearn.decomposition


def plot_explained_variance(
    pca: sklearn.decomposition.PCA, ax: Optional[matplotlib.axes.Axes] = None
) -> matplotlib.axes.Axes:
    """Plot explained variance of PCA from scikit-learn."""
    exp_var = pd.Series(pca.explained_variance_ratio_).to_frame("explained variance")
    exp_var.index += 1  # start at 1
    exp_var["explained variance (cummulated)"] = exp_var["explained variance"].cumsum()
    exp_var.index.name = "PC"
    ax = exp_var.plot(ax=ax)
    return ax
