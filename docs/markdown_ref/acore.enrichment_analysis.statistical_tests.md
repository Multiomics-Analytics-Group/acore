# acore.enrichment_analysis.statistical_tests namespace

## Submodules

## acore.enrichment_analysis.statistical_tests.fisher module

Run fisher’s exact test on two groups using [scipy.stats.fisher_exact](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html).

### run_fisher(group1: list[int], group2: list[int], alternative: str = 'two-sided') → tuple[float, float]

Run fisher’s exact test on two groups using [scipy.stats.fisher_exact](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html).

Example:

```default
# annotated   not-annotated
# group1      a               b
# group2      c               d


odds, pvalue = stats.fisher_exact(group1=[a, b],
                                  group2 =[c, d]
                )
```

## acore.enrichment_analysis.statistical_tests.kolmogorov_smirnov module

### run_kolmogorov_smirnov(dist1: list[float], dist2: list[float], alternative: str = 'two-sided') → tuple[float, float]

Compute the Kolmogorov-Smirnov statistic on 2 samples.
See [scipy.stats.ks_2samp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html)

* **Parameters:**
  * **dist1** (*list*) – sequence of 1-D ndarray (first distribution to compare)
    drawn from a continuous distribution
  * **dist2** (*list*) – sequence of 1-D ndarray (second distribution to compare)
    drawn from a continuous distribution
  * **alternative** (*str*) – defines the alternative hypothesis (default is ‘two-sided’):
    \* **‘two-sided’**
    \* **‘less’**
    \* **‘greater’**
* **Returns:**
  statistic float and KS statistic pvalue float Two-tailed p-value.

Example:

```default
result = run_kolmogorov_smirnov(dist1, dist2, alternative='two-sided')
```
