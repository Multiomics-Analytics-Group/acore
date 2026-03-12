# acore.enrichment_analysis.statistical_tests namespace

## Submodules

## acore.enrichment_analysis.statistical_tests.fisher module

Run fisher’s exact test on two groups using [scipy.stats.fisher_exact](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html).

### run_fisher(group1: [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)], group2: [list](https://docs.python.org/3/library/stdtypes.html#list)[[int](https://docs.python.org/3/library/functions.html#int)], alternative: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'two-sided') → [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)[[float](https://docs.python.org/3/library/functions.html#float), [float](https://docs.python.org/3/library/functions.html#float)]

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

### run_kolmogorov_smirnov(dist1: [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)], dist2: [list](https://docs.python.org/3/library/stdtypes.html#list)[[float](https://docs.python.org/3/library/functions.html#float)], alternative: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'two-sided') → [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)[[float](https://docs.python.org/3/library/functions.html#float), [float](https://docs.python.org/3/library/functions.html#float)]

Compute the Kolmogorov-Smirnov statistic on 2 samples.
See [scipy.stats.ks_2samp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html)

* **Parameters:**
  * **dist1** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – sequence of 1-D ndarray (first distribution to compare)
    drawn from a continuous distribution
  * **dist2** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – sequence of 1-D ndarray (second distribution to compare)
    drawn from a continuous distribution
  * **alternative** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – defines the alternative hypothesis (default is ‘two-sided’):
    \* **‘two-sided’**
    \* **‘less’**
    \* **‘greater’**
* **Returns:**
  statistic float and KS statistic pvalue float Two-tailed p-value.

Example:

```default
result = run_kolmogorov_smirnov(dist1, dist2, alternative='two-sided')
```
