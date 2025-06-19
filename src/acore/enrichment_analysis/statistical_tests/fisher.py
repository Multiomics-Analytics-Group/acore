"""Run fisher's exact test on two groups using `scipy.stats.fisher_exact`_."""

from __future__ import annotations

from scipy import stats


def run_fisher(
    group1: list[int],
    group2: list[int],
    alternative: str = "two-sided",
) -> tuple[float, float]:
    """Run fisher's exact test on two groups using `scipy.stats.fisher_exact`_.

    .. _scipy.stats.fisher_exact: https://docs.scipy.org/doc/scipy/reference/generated/\
scipy.stats.fisher_exact.html

    Example::

        # annotated   not-annotated
        # group1      a               b
        # group2      c               d


        odds, pvalue = stats.fisher_exact(group1=[a, b],
                                          group2 =[c, d]
                        )
    """

    odds, pvalue = stats.fisher_exact([group1, group2], alternative)

    return (odds, pvalue)
