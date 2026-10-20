import numpy as np
import pandas as pd

from acore.enrichment_analysis.statistical_tests.fisher import run_fisher
from acore.multiple_testing import compute_efdr


def _compute_efdr_for_enrichment(
    data: pd.DataFrame,
    observed_pvalues: list,
    observed_terms: list,
    foreground_id: str,
    background_id: str,
    foreground_pop: int,
    background_pop: int,
    min_detected_in_set: int,
    annotation_col: str,
    group_col: str,
    identifier_col: str,
    permutations: int,
    alpha: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Internal helper that generates permuted p-values and calls :func:`compute_efdr`.

    For each permutation the group labels are randomly shuffled across the unique
    annotated features (preserving the original group sizes), so that a feature
    keeps one label across all of its annotation terms. Fisher's test is re-run for every
    annotation term present in the observed result, and p-values are collected.
    Terms that do not meet ``min_detected_in_set`` in a given permutation are
    assigned p-value 1.0.
    """
    # one label per feature, as a feature can be annotated with multiple terms
    feature_groups = data.drop_duplicates(identifier_col).set_index(identifier_col)[
        group_col
    ]
    rng = np.random.default_rng()

    permuted_pvalues = []
    for _ in range(permutations):
        perm_feature_groups = pd.Series(
            rng.permutation(feature_groups.to_numpy()), index=feature_groups.index
        )
        perm_data = data.copy()
        perm_data[group_col] = perm_data[identifier_col].map(perm_feature_groups)

        perm_countsdf = (
            perm_data.groupby([annotation_col, group_col])
            .agg(["count"])[(identifier_col, "count")]
            .reset_index()
        )
        perm_countsdf.columns = [annotation_col, group_col, "count"]

        perm_pvals = []
        for term in observed_terms:
            perm_counts = perm_countsdf[perm_countsdf[annotation_col] == term]
            num_fg = perm_counts.loc[
                perm_counts[group_col] == foreground_id, "count"
            ].sum()
            num_bg = perm_counts.loc[
                perm_counts[group_col] == background_id, "count"
            ].sum()
            if num_fg >= min_detected_in_set:
                _, pval = run_fisher(
                    [num_fg, foreground_pop - num_fg],
                    [num_bg, background_pop - foreground_pop - num_bg],
                )
            else:
                pval = 1.0
            perm_pvals.append(pval)

        permuted_pvalues.append(np.array(perm_pvals))

    return compute_efdr(observed_pvalues, permuted_pvalues, alpha=alpha)
