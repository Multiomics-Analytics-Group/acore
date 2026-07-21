# acore.permutation_test package

### paired_permutation(cond1: [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray), cond2: [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray), metric: [str](https://docs.python.org/3/library/stdtypes.html#str) = 't-statistic', n_permutations: [int](https://docs.python.org/3/library/functions.html#int) = 10000, rng: [Generator](https://numpy.org/doc/stable/reference/random/generator.html#numpy.random.Generator) = Generator(PCG64) at 0x7F064A9A7220, \*\*kwargs) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)

Perform a permutation test for paired samples.

* **Parameters:**
  * **cond1** (*np.ndarray*) – First condition (paired samples).
  * **cond2** (*np.ndarray*) – Second condition (paired samples).
  * **metric** ([*str*](https://docs.python.org/3/library/stdtypes.html#str) *or* *callable* *,* *optional*) – Metric to compute (‘t-statistic’, ‘mean’, ‘median’, or a custom function).
  * **n_permutations** ([*int*](https://docs.python.org/3/library/functions.html#int) *,* *optional*) – Number of permutations to perform (default is 10000).
  * **rng** (*np.random.Generator* *,* *optional*) – Random number generator (default is np.random.default_rng(seed=12345)).
  * **\*\*kwargs** – Additional arguments passed to the metric function.
* **Returns:**
  Dictionary with keys:
  - ‘metric’: Metric function used.
  - ‘observed’: Observed metric value.
  - ‘p_value’: Permutation test p-value (np.nan if degenerate).
* **Return type:**
  [dict](https://docs.python.org/3/library/stdtypes.html#dict)

### chi2_permutation(\*groups, n_permutations: [int](https://docs.python.org/3/library/functions.html#int) = 10000, rng: [Generator](https://numpy.org/doc/stable/reference/random/generator.html#numpy.random.Generator) = Generator(PCG64) at 0x7F064A9A6B20) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)

Perform a permutation test for categorical data using the chi-squared statistic.

* **Parameters:**
  * **\*groups** (*array-like*) – Arrays representing categorical groups.
  * **n_permutations** ([*int*](https://docs.python.org/3/library/functions.html#int) *,* *optional*) – Number of permutations to perform (default is 10000).
  * **rng** (*np.random.Generator* *,* *optional*) – Random number generator (default is np.random.default_rng(seed=12345)).
* **Returns:**
  Dictionary with keys:
  - ‘observed’: Observed chi-squared test result.
  - ‘p_value’: Permutation test p-value.
* **Return type:**
  [dict](https://docs.python.org/3/library/stdtypes.html#dict)

### indep_permutation(group1: [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray), group2: [ndarray](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray), metric: [str](https://docs.python.org/3/library/stdtypes.html#str) = 't-statistic', n_permutations: [int](https://docs.python.org/3/library/functions.html#int) = 10000, rng: [Generator](https://numpy.org/doc/stable/reference/random/generator.html#numpy.random.Generator) = Generator(PCG64) at 0x7F064A9A73E0, \*\*kwargs) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)

Perform a permutation test for independent samples.

* **Parameters:**
  * **group1** (*np.ndarray*) – First group of samples.
  * **group2** (*np.ndarray*) – Second group of samples.
  * **metric** ([*str*](https://docs.python.org/3/library/stdtypes.html#str) *or* *callable* *,* *optional*) – Metric to compute (‘t-statistic’, ‘anova’, ‘mean’, ‘median’, or a custom function).
  * **n_permutations** ([*int*](https://docs.python.org/3/library/functions.html#int) *,* *optional*) – Number of permutations to perform (default is 10000).
  * **rng** (*np.random.Generator* *,* *optional*) – Random number generator (default is np.random.default_rng(seed=12345)).
  * **\*\*kwargs** – Additional arguments passed to the metric function.
* **Returns:**
  Dictionary with keys:
  - ‘metric’: Metric function used.
  - ‘observed’: Observed metric value.
  - ‘p_value’: Permutation test p-value.
* **Return type:**
  [dict](https://docs.python.org/3/library/stdtypes.html#dict)

## Submodules

## acore.permutation_test.internal_functions module

## acore.permutation_test.jaccard module

### jaccard_similarity(set1: [set](https://docs.python.org/3/library/stdtypes.html#set), set2: [set](https://docs.python.org/3/library/stdtypes.html#set)) → [float](https://docs.python.org/3/library/functions.html#float)

Compute the Jaccard similarity between two sets.

* **Parameters:**
  * **set1** ([*set*](https://docs.python.org/3/library/stdtypes.html#set)) – First set of nodes.
  * **set2** ([*set*](https://docs.python.org/3/library/stdtypes.html#set)) – Second set of nodes.
* **Returns:**
  Jaccard similarity coefficient,
  which is the size of the intersection divided by
  the size of the union of the two sets.
* **Return type:**
  [float](https://docs.python.org/3/library/functions.html#float)

### Example

```pycon
>>> set1 = {1, 2, 3}
>>> set2 = {2, 3, 4}
>>> jaccard_similarity(set1, set2)
0.5
```

### avg_jaccard(group: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)) → [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)

> Compute the average pairwise Jaccard similarity within a group of graphs.

> group
> : An iterator of iterable objects, such as a list of sets.

Returns
: tuple
  : A tuple containing the average Jaccard similarity and its standard deviation.
  <br/>
  ```pycon
  >>> group = [{1, 2, 3}, {2, 3, 4}, {3, 4, 5}]
  >>> avg_jaccard(group)
  (0.3333333333333333, 0.0)
  ```

### btwn_jaccard(group1: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable), group2: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)) → [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)

Compute the average Jaccard similarity between two groups of graphs.

* **Parameters:**
  * **group1** (*Iterable*) – An iterator of iterable objects, such as a list of sets.
  * **group2** (*Iterable*) – Another iterator of iterable objects.
* **Returns:**
  A tuple containing the average Jaccard similarity and its standard deviation.
* **Return type:**
  [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)

### Example

```pycon
>>> group1 = [{1, 2, 3}, {2, 3, 4}]
>>> group2 = [{3, 4, 5}, {4, 5, 6}]
>>> btwn_jaccard(group1, group2)
(0.3333333333333333, 0.0)
```
