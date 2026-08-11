"""Filtering annotations by keywords.

from acore.io.uniprot import fetch_annotations, process_annotations, filter_annotations
fields = "accession,go_p,go_c,go_f"
uniprot_ids = ["P12345", "Q67890", "P05067", "A1B2C3"]
df = fetch_annotations(uniprot_ids, fields=fields)
df = process_annotations(df, fields=fields)
# first keyword does not exist
filtered = filter_annotations(df, keywords=["apoptosis", "mitochondr", "synapse"])
"""

from __future__ import annotations

import re as _re

import pandas as pd

re_IGNORECASE = _re.IGNORECASE

__all__ = ["filter_annotations"]


def filter_annotations(
    annotations: pd.DataFrame,
    keywords: list[str],
    case_sensitive: bool = False,
) -> pd.DataFrame:
    """Filter a long-format annotations DataFrame to rows whose ``annotation``
    column contains **any** of the given keywords.

    Parameters
    ----------
    annotations : pd.DataFrame
        Long-format DataFrame with at least an ``annotation`` column (as
        returned by :func:`process_annotations`).
    keywords : list[str]
        Keywords to search for.  A row is kept when at least one keyword
        appears in the annotation string.  Empty strings are ignored.
    case_sensitive : bool, optional
        Whether the search is case-sensitive.  Default ``False``.

    Returns
    -------
    pd.DataFrame
        Filtered copy of *annotations* where at least one keyword matched.
        Returns an empty DataFrame with the same columns when no keywords
        are provided or none match.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     "identifier": ["P1", "P1", "P2"],
    ...     "source": ["go_p", "go_p", "go_p"],
    ...     "annotation": [
    ...         "cell apoptosis [GO:0042981]",
    ...         "cell division",
    ...         "mitochondria",
    ...     ],
    ... })
    >>> result = filter_annotations(df, keywords=["apoptosis", "mitochondr"])
    >>> list(result["annotation"])
    ['cell apoptosis [GO:0042981]', 'mitochondria']
    >>> filter_annotations(df, keywords=[]).shape[0]
    0
    """
    # Drop empty / whitespace-only keywords
    clean_keywords = [kw.strip() for kw in keywords if kw.strip()]
    if not clean_keywords:
        return annotations.iloc[0:0].copy()

    flags = 0 if case_sensitive else re_IGNORECASE
    pattern = "|".join(_escape(kw) for kw in clean_keywords)
    mask = annotations["annotation"].str.contains(pattern, flags=flags, na=False)
    return annotations[mask].copy()


def _escape(keyword: str) -> str:
    return _re.escape(keyword)
