# acore.multiple_testing package

### apply_pvalue_correction(pvalues: [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray), alpha: [float](https://docs.python.org/3/library/functions.html#float) = 0.05, method: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'bonferroni') → [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)[[ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray), [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)]

Performs p-value correction using the specified method as in
[statsmodels.stats.multitest.multipletests](https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html).

* **Parameters:**
  * **pvalues** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – et of p-values of the individual tests.
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate.
  * **method** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – 

    method of p-value correction:
    - ’bonferroni’ : one-step correction
    - ’sidak’ : one-step correction
    - ’holm-sidak’ : step down method using Sidak adjustments
    - ’holm’ : step-down method using Bonferroni adjustments
    - ’simes-hochberg’ : step-up method (independent)
    - ’hommel’ : closed method based on Simes tests (non-negative)
    - ’fdr_bh’ : Benjamini/Hochberg (non-negative)
    - ’fdr_by’ : Benjamini/Yekutieli (negative)
    - ’fdr_tsbh’ : two stage fdr correction (non-negative)
    - ’fdr_tsbky’ : two stage fdr correction (non-negative)
* **Returns:**
  Tuple with two numpy.array\`s, boolen for rejecting H0 hypothesis
  and float for adjusted p-value. Can contain missing values if \`pvalues
  contain missing values.

Example:

```default
result = apply_pvalue_correction(pvalues, alpha=0.05, method='bonferroni')
```

### apply_pvalue_fdrcorrection(pvalues, alpha=0.05, method='indep')

Performs p-value correction for false discovery rate. For more information visit [https://www.statsmodels.org/devel/generated/statsmodels.stats.multitest.fdrcorrection.html](https://www.statsmodels.org/devel/generated/statsmodels.stats.multitest.fdrcorrection.html).

* **Parameters:**
  * **pvalues** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – et of p-values of the individual tests.
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate.
  * **method** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method of p-value correction (‘indep’, ‘negcorr’).
* **Returns:**
  Tuple with two arrays, boolen for rejecting H0 hypothesis and float for adjusted p-value.

Example:

```default
result = apply_pvalue_fdrcorrection(pvalues, alpha=0.05, method='indep')
```

### apply_pvalue_twostage_fdrcorrection(pvalues, alpha=0.05, method='bh')

Iterated two stage linear step-up procedure with estimation of number of true hypotheses. For more information visit [https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.fdrcorrection_twostage.html](https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.fdrcorrection_twostage.html).

* **Parameters:**
  * **pvalues** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – et of p-values of the individual tests.
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate.
  * **method** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – method of p-value correction (‘bky’, ‘bh’).
* **Returns:**
  Tuple with two arrays, boolen for rejecting H0 hypothesis and float for adjusted p-value.

Example:

```default
result = apply_pvalue_twostage_fdrcorrection(pvalues, alpha=0.05, method='bh')
```

### apply_pvalue_permutation_fdrcorrection(df, observed_pvalues, group, alpha=0.05, permutations=50)

This function applies multiple hypothesis testing correction using a permutation-based false discovery rate approach.

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and features as columns.
  * **oberved_pvalues** – pandas Series with p-values calculated on the originally measured data.
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the column containing group identifiers.
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate. Values velow alpha are considered significant.
  * **permutations** ([*int*](https://docs.python.org/3/library/functions.html#int)) – number of permutations to be applied.
* **Returns:**
  Pandas dataframe with adjusted p-values and rejected columns.

Example:

```default
result = apply_pvalue_permutation_fdrcorrection(df, observed_pvalues, group='group', alpha=0.05, permutations=50)
```

### calculate_anova(df, column, group='group')

Calculates one-way ANOVA using pingouin.

* **Parameters:**
  * **df** – pandas dataframe with group as rows and protein identifier as column
  * **column** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the column in df to run ANOVA on
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with group identifiers
* **Returns:**
  Tuple with t-statistics and p-value.

### get_counts_permutation_fdr(value, random, observed, n, alpha)

Calculates local FDR values (q-values) by computing the fraction of accepted hits from the permuted data over accepted hits from the measured data normalized by the total number of permutations.

* **Parameters:**
  * **value** ([*float*](https://docs.python.org/3/library/functions.html#float)) – computed p-value on measured data for a feature.
  * **random** ([*numpy.ndarray*](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)) – p-values computed on the permuted data.
  * **observed** – pandas Series with p-values calculated on the originally measured data.
  * **n** ([*int*](https://docs.python.org/3/library/functions.html#int)) – number of permutations to be applied.
  * **alpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – error rate. Values velow alpha are considered significant.
* **Returns:**
  Tuple with q-value and boolean for H0 rejected.

Example:

```default
result = get_counts_permutation_fdr(value, random, observed, n=250, alpha=0.05)
```

### get_max_permutations(df, group='group')

Get maximum number of permutations according to number of samples.

* **Parameters:**
  * **df** – pandas dataframe with samples as rows and protein identifiers as columns
  * **group** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column with group identifiers
* **Returns:**
  Maximum number of permutations.
* **Return type:**
  [int](https://docs.python.org/3/library/functions.html#int)

### correct_pairwise_ttest(df, alpha, correction='fdr_bh')
