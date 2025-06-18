"""Put unique features into foreground, background or assign nan."""

from __future__ import annotations

import numpy as np
import pandas as pd


def annotate_features(
    features: pd.Series,
    in_foreground: set[str] | list[str],
    in_background: set[str] | list[str],
) -> pd.Series:
    """
    Annotate features as foreground or background based on their presence in the
    foreground and background lists.

    :param features: pandas.Series with features and their annotations.
    :param in_foreground: list of features identifiers in the foreground.
    :type in_foreground: set or list-like
    :param in_background: list of features identifiers in the background.
    :type in_background: set or list-like
    :return: pandas.Series containing 'foreground' or 'background'.
             missing values are preserved.

    Example::

        result = _annotate_features(features, in_foreground, in_background)
    """
    in_either_or = features.isin(in_foreground) | features.isin(in_background)
    res = (
        features.where(in_either_or, np.nan)
        .mask(features.isin(in_foreground), "foreground")
        .mask(features.isin(in_background), "background")
    )
    return res
