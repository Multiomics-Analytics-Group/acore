"""Module for filtering metabolomics feature table."""

import logging

import pandas as pd
import scipy

logger = logging.getLogger(__name__)


def filter_by_missingness(
    data: pd.DataFrame,
    percent: int = 80,
    method: str = "classic",
    samples: list | None = None,
    groups: dict | None = None,
):
    """
    Implementation of the 80%-rule.

    If there are more than 20% of values (intensities) missing for one feature,
    this feature will get removed.

    :param data: pandas data frame with samples as rows and features as columns.
    :param percent: percentage chosen for filtering. The default is 80%, meaning that
        at least 80% of the values of every feature need to be present in order for this
        feature to be retained.
    :param method: str that is either "classic" or "modified".
        If "classic", all samples are considered for each feature. Samples are taken from the
        "samples" parameter and should not include controls or QCs.
        If "modified", conditions are separated when calculating the percentage of
        missingness. A feature is retained if at least `percent`% of its values are
        present in ANY one condition. This allows condition-specific features (e.g.
        present in treatment but missing in control) to be retained.
    :param samples: list of row index labels (from data.index) identifying the biological
        sample rows, e.g. ["S1", "S2", "S3"]. Required when method="classic". Should not
        include control or QC samples.
    :param groups: dict mapping condition name to a list of row index labels belonging to
        that condition, e.g. {"treatment": ["S1", "S2", "S3"], "control": ["S4", "S5", "S6"]}.
        Required when method="modified", ignored otherwise. QCs and blanks are excluded
        by simply not including them in the dict.
    """
    if not 0 < percent < 100:
        logger.warning(
            f"The percentage value must be between 0 and 100. User input: {percent}"
        )

    df = data.copy()

    threshold = percent / 100

    if method == "classic":
        if samples is None:
            raise ValueError("sample names must be provided when method='classic'")
        keep = df.loc[samples].notna().mean() >= threshold
    elif method == "modified":
        if groups is None:
            raise ValueError("groups must be provided when method='modified'")
        keep = pd.Series(False, index=df.columns)
        for condition_samples in groups.values():
            group_data = df.loc[condition_samples]
            keep |= group_data.notna().mean() >= threshold
    else:
        raise ValueError(f"method must be 'classic' or 'modified', got '{method}'")

    return df.loc[:, keep]


def filter_cv(
    data: pd.DataFrame,
    samples: list,
    qcs: list,
):
    """
    Implementation of coefficient of variation (CV)-based filtering.

    Features are removed when their CV across biological samples is smaller than their CV
    across QC samples, meaning analytical noise exceeds biological variability.

    :param data: pandas data frame with samples as rows and features as columns.
    :param samples: list of row index labels (from data.index) identifying the
        biological sample rows, e.g. ["S1", "S2", "S3"].
    :param qcs: list of row index labels identifying the quality control rows,
        e.g. ["QC1", "QC2", "QC3"].
    """

    def cv(df):
        # scipy.stats.variation(df, ddof=1)
        return df.std() / df.mean()

    if len(qcs) < 2:
        raise ValueError(
            f"You need more than 1 QC sample to apply this filtering method. Got {len(qcs)}."
        )

    df = data.copy()
    cv_qcs = cv(df.loc[qcs])

    undefined = cv_qcs.isna() | cv_qcs.isin([float("inf"), float("-inf")])
    if undefined.any():
        logger.warning(
            "CV is undefined for %d feature(s) in QC samples (zero or near-zero mean): "
            "%s. These features will be dropped. Consider running filter_by_missingness first.",
            undefined.sum(),
            list(df.columns[undefined]),
        )

    cv_samples = cv(df.loc[samples])
    keep = cv_samples >= cv_qcs

    return df.loc[:, keep]


def filter_blanks(
    data: pd.DataFrame,
    blanks: list,
    samples: list,
    threshold: float = 0.5,
):
    """
    Filtering out features that show up in the blanks control.

    The mean intensity scores are calculated per-feature within the
    blanks and the samples. If the ratio of a feature's mean intensity in the blanks
    to its mean intensity in the samples is more than half (per default), the feature
    gets removed. It is assumed to have potentially contaminated the instrument, so
    the measurements in the samples cannot be trusted to be biologically relevant.

    :param data: pandas DataFrame containing data with samples as rows and features as columns
    :param blanks: list of row index labels (from data.index) identifying the blanks
        measurement rows, e.g. ["Blank1", "Blank2"]
    :param samples: list of row index labels (from data.index) identifying the biological
        sample rows, e.g. ["S1", "S2", "S3"]
    :param threshold: optional ratio used as a threshold to determine whether the detected
        intensities in blanks are too high in comparison with sample intensities.
        Defaults to 0.5, but can be adjusted based on data and stringency.
    """
    df = data.copy()

    mean_blank = df.loc[blanks].mean(axis=0)
    mean_sample = df.loc[samples].mean(axis=0)

    epsilon = 1e-9
    ratio = (mean_blank + epsilon) / (mean_sample + epsilon)

    features_to_remove = df.columns[ratio > threshold]

    return df.drop(columns=features_to_remove)
