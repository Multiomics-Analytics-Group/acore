"""Plot metrics for binary classification."""

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from njab.sklearn.types import Results, ResultsSplit

LIMITS = (-0.05, 1.05)


def plot_split_auc(
    result: ResultsSplit, name: str, ax: matplotlib.axes.Axes
) -> matplotlib.axes.Axes:
    """Add receiver operation curve to ax of a split of the data."""
    col_name = f"{name} (auc: {result.auc:.3f})"
    roc = pd.DataFrame(result.roc, index="fpr tpr cutoffs".split()).rename(
        {"tpr": col_name}
    )
    ax = roc.T.plot(
        "fpr",
        col_name,
        xlabel="false positive rate",
        ylabel="true positive rate",
        style=".-",
        ylim=LIMITS,
        xlim=LIMITS,
        ax=ax,
    )
    return ax


# ! should be roc
def plot_auc(
    results: Results,
    ax: matplotlib.axes.Axes = None,
    label_train="train",
    label_test="test",
    **kwargs,
) -> matplotlib.axes.Axes:
    """Plot ROC curve for train and test data."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, **kwargs)
    ax = plot_split_auc(results.train, f"{label_train}", ax)
    ax = plot_split_auc(results.test, f"{label_test}", ax)
    return ax


def plot_split_prc(
    result: ResultsSplit, name: str, ax: matplotlib.axes.Axes
) -> matplotlib.axes.Axes:
    """Add precision recall curve to ax of a split of the data."""
    col_name = f"{name} (aps: {result.aps:.3f})"
    roc = pd.DataFrame(result.prc, index="precision recall cutoffs".split()).rename(
        {"precision": col_name}
    )
    ax = roc.T.plot(
        "recall",
        col_name,
        xlabel="true positive rate",
        ylabel="precision",
        style=".-",
        ylim=LIMITS,
        xlim=LIMITS,
        ax=ax,
    )
    return ax


def plot_prc(
    results: ResultsSplit,
    ax: matplotlib.axes.Axes = None,
    label_train="train",
    label_test="test",
    **kwargs,
):
    """Plot precision recall curve for train and test data."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, **kwargs)
    ax = plot_split_prc(results.train, f"{label_train}", ax)
    ax = plot_split_prc(results.test, f"{label_test}", ax)
    return ax
