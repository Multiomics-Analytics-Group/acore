import numpy as np
import pandas as pd
from scipy import sparse

from acore.enrichment_analysis.statistical_tests.fisher import run_fisher
from acore.multiple_testing import compute_efdr


def _compute_efdr_for_enrichment(
    data: pd.DataFrame,
    observed_pvalues: list,
    foreground_pop: int,
    background_pop: int,
    min_detected_in_set: int,
    annotation_col: str,
    identifier_col: str,
    permutations: int,
    alpha: float,
    alternative: str = "two-sided",
    seed: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Internal helper that generates permuted p-values and calls :func:`compute_efdr`.

    In each permutation a random foreground set of the original size
    (``foreground_pop``) is drawn from *all* ``background_pop`` features, i.e. also
    from features without any annotation. The permuted foreground is therefore
    consistent with the null hypothesis of the (hypergeometric) Fisher's test and
    does not condition on the number of annotated features in the foreground. A
    feature keeps one label across all of its annotation terms. This follows the
    eFDR procedure of the mulea R package (Turek et al. 2024).

    Fisher's test is re-run for *every* annotation term in ``data`` (not only for
    the terms tested in the observed data) and all permuted p-values are pooled to
    estimate the number of false positives expected at each p-value threshold. Terms
    that do not meet ``min_detected_in_set`` in a given permutation are assigned
    p-value 1.0.

    :param data: annotation of the features in the foreground or background
        (one row per feature and annotation term).
    :param observed_pvalues: p-values of the tested terms in the observed data.
    :param seed: seed of the random number generator for reproducible permutations.
    :param alternative: alternative hypothesis of Fisher's test, see
        :func:`acore.enrichment_analysis.run_fisher`.
    """
    data = data[[identifier_col, annotation_col]].drop_duplicates()
    feature_codes, _ = pd.factorize(data[identifier_col])
    term_codes, terms = pd.factorize(data[annotation_col])
    n_features = feature_codes.max() + 1
    # features x terms membership matrix
    membership = sparse.csr_matrix(
        (np.ones(len(data), dtype=np.int32), (feature_codes, term_codes)),
        shape=(n_features, len(terms)),
    )
    # number of annotated features per term (foreground and background)
    n_in_term = np.asarray(membership.sum(axis=0)).ravel()

    rng = np.random.default_rng(seed)
    cache: dict[int, float] = {}
    stride = foreground_pop + 1

    def _pvalue(key: int) -> float:
        if key not in cache:
            n_in, num_fg = divmod(key, stride)
            _, cache[key] = run_fisher(
                [num_fg, foreground_pop - num_fg],
                [n_in - num_fg, background_pop - foreground_pop - (n_in - num_fg)],
                alternative=alternative,
            )
        return cache[key]

    permuted_pvalues = []
    for _ in range(permutations):
        # annotated features are slots 0..n_features-1 of the background population
        drawn = rng.choice(background_pop, size=foreground_pop, replace=False)
        in_foreground = np.zeros(n_features, dtype=np.int32)
        in_foreground[drawn[drawn < n_features]] = 1
        num_fg = membership.T.dot(in_foreground)
        keys, inverse = np.unique(n_in_term * stride + num_fg, return_inverse=True)
        pvalues = np.array([_pvalue(int(key)) for key in keys])[inverse]
        pvalues[num_fg < min_detected_in_set] = 1.0
        permuted_pvalues.append(pvalues)

    return compute_efdr(observed_pvalues, permuted_pvalues, alpha=alpha)
