"""Formatting related functions for test tables."""

import numpy as np
import pandas as pd

from acore.multiple_testing import (
    apply_pvalue_correction,
    apply_pvalue_permutation_fdrcorrection,
    get_max_permutations,
)


def format_anova_table(
    df: pd.DataFrame,
    aov_results: list[tuple],
    pairwise_results,
    pairwise_cols,
    group: str,
    permutations: int,
    alpha: float,
    correction: str,
) -> pd.DataFrame:
    """
    Performs p-value correction (permutation-based and FDR) and converts pandas dataframe
    into final format.

    :param df: pandas dataframe with samples as rows and protein identifiers as columns
               (with additional columns 'group', 'sample' and 'subject').
    :param list[tuple] aov_results: list of tuples with anova results (one tuple per feature).
    :param list[dataframes] pairwise_results: list of pandas dataframes with
                                              posthoc tests results
    :param str group: column with group identifiers
    :param float alpha: error rate for multiple hypothesis correction
    :param int permutations: number of permutations used to estimate false discovery rates
    :return: Pandas dataframe
    """
    columns = ["identifier", "dfk", "dfn", "F-statistics", "pvalue"]
    scores = pd.DataFrame(aov_results, columns=columns)
    scores = scores.set_index("identifier")
    corrected = False
    # FDR correction
    if permutations > 0:
        max_perm = get_max_permutations(df, group=group)
        if max_perm >= 10:
            permutations = min(permutations, max_perm)
            observed_pvalues = scores.pvalue
            count = apply_pvalue_permutation_fdrcorrection(
                df,
                observed_pvalues,
                group=group,
                alpha=alpha,
                permutations=permutations,
            )
            scores = scores.join(count)
            scores["correction"] = f"permutation FDR ({permutations} perm)"
            corrected = True

    if not corrected:
        _, padj = apply_pvalue_correction(
            scores["pvalue"].tolist(), alpha=alpha, method=correction
        )
        scores["correction"] = correction
        scores["padj"] = padj
        corrected = True
    # ! scores is omnibus result
    res = pd.DataFrame(pairwise_results, columns=pairwise_cols).set_index("identifier")
    if not res.empty:
        res = res.join(scores[["F-statistics", "pvalue", "padj"]].astype("float"))
        res["correction"] = scores["correction"]
    else:
        res = scores
        res["log2FC"] = np.nan

    res = res.reset_index()
    if {"group1", "group2"}.issubset(res.columns):
        res["posthoc comparison"] = (
            res["group1"].astype(str) + "~~" + res["group2"].astype(str)
        )
    res["rejected"] = res["padj"] < alpha
    # ! double check this
    if "posthoc pvalue" in res.columns:
        res["-log10 posthoc pvalue"] = [
            -np.log10(x) for x in res["posthoc pvalue"].values
        ]
    if "pvalue" in res.columns:
        res["-log10 pvalue"] = [-np.log10(x) for x in res["pvalue"].values]

    return res
