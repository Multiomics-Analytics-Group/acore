# %% [markdown]
# (nb_ref_efdr_mulea)=
# # Empirical FDR for overrepresentation analysis: comparison with mulea
# The [mulea](https://github.com/ELTEbioinformatics/mulea) R package
# ([Turek et al. 2024](https://doi.org/10.1093/bib/bbae301)) introduced an empirical false
# discovery rate (eFDR) for overrepresentation analysis (ORA) which accounts for the
# overlap between annotations (here: transcription factor regulons). Setting
# `correction="efdr"` in the enrichment functions of `acore` uses the same procedure:
#
# 1. test all annotations with a (one-sided) Fisher's exact test, i.e. a hypergeometric
#    test
# 2. repeat the test for all annotations on `permutations` random foregrounds of the same
#    size drawn from the background
# 3. at the p-value $t$ of an annotation, the eFDR is the average number of annotations
#    with a p-value $\leq t$ in the permutations, divided by the number of annotations
#    with a p-value $\leq t$ in the observed data
#
# Here we reproduce the example from the mulea README: differentially expressed genes of
# *E. coli* treated with ciprofloxacin
# ([GSE55662](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55662)) and the
# transcription factor regulons of [RegulonDB](https://regulondb.ccg.unam.mx/) as
# annotation.

# %% tags=["hide-output"]
# %pip install acore

# %%
from urllib.request import urlopen

import pandas as pd

import acore.enrichment_analysis
from acore.multiple_testing import apply_pvalue_correction

# %% [markdown]
# # Load the data used in the mulea README
# - `target_set`: significantly regulated genes (foreground)
# - `background_set`: all genes investigated in the experiment
# - the transcription factor ontology in the GMT format

# %%
MULEA = "https://raw.githubusercontent.com/ELTEbioinformatics/mulea/master/inst/extdata"
GMT = (
    "https://raw.githubusercontent.com/ELTEbioinformatics/GMT_files_for_mulea/main/"
    "GMT_files/Escherichia_coli_83333/"
    "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"
)


def read_lines(url: str) -> list[str]:
    with urlopen(url) as response:
        return response.read().decode().splitlines()


target_set = set(read_lines(f"{MULEA}/target_set.txt"))
background_set = read_lines(f"{MULEA}/background_set.txt")
print(
    f"{len(target_set)} genes in the target set, {len(background_set)} in the background"
)

# %% [markdown]
# Parse the GMT file (one transcription factor per line: id, name, target genes) and
# keep the same terms as in the mulea example. Note that `filter_ontology` of mulea
# uses strict inequalities, so `min_nr_of_elements = 3` and `max_nr_of_elements = 400`
# keep terms with 4 to 399 genes.

# %%
annotation = []
for line in read_lines(GMT):
    if line.startswith("#") or not line.strip():
        continue
    term_id, _name, *genes = line.split("\t")
    genes = [gene for gene in genes if gene]
    if 3 < len(genes) < 400:
        annotation.extend((term_id, gene) for gene in set(genes))
annotation = pd.DataFrame(annotation, columns=["annotation", "identifier"])
print(f"{annotation['annotation'].nunique()} transcription factors in the ontology")
annotation

# %% [markdown]
# For `acore` the data needs to be in the form of a regulation table with one row per
# feature of the background and a column `rejected`, which is `True` for the target set.

# %%
regulation_data = pd.DataFrame({"identifier": background_set})
regulation_data["rejected"] = regulation_data["identifier"].isin(target_set)
regulation_data["rejected"].value_counts()

# %% [markdown]
# # Overrepresentation analysis with eFDR
# Enrichment analysis uses a one-sided test for overrepresentation, as mulea, with
# `min_detected_in_set=1`, and the 10000 permutations recommended by mulea.
# Fixing the `seed` makes the result reproducible.

# %%
efdr = acore.enrichment_analysis.run_regulation_enrichment(
    regulation_data=regulation_data,
    annotation=annotation.copy(),
    min_detected_in_set=1,
    correction="efdr",
    permutations=10_000,
    seed=1,
).set_index("terms")
efdr["rejected"] = efdr["rejected"].astype(bool)
efdr.query("rejected")[["foreground", "pvalue", "padj"]].round(5)

# %% [markdown]
# # Comparison with the result reported by mulea
# The table below is the result of the mulea README (`ora_results` filtered for
# `eFDR < 0.05`, 10000 permutations, `random_seed = 1`).

# %% tags=["hide-input"]
mulea_reference = pd.DataFrame(
    {
        "terms": [
            "FNR",
            "LexA",
            "SoxS",
            "Rob",
            "DnaA",
            "FadR",
            "NsrR",
            "ArcA",
            "IHF",
            "MarA",
        ],
        "foreground": [26, 14, 7, 5, 4, 5, 8, 12, 14, 5],
        "pvalue": [
            3e-7,
            0.0,
            0.0001615,
            0.0004717,
            0.0006281,
            0.0003692,
            0.0010478,
            0.0032001,
            0.0070758,
            0.0066068,
        ],
        "eFDR": [
            0.0,
            0.0,
            0.0036667,
            0.00512,
            0.0052,
            0.0056,
            0.0073714,
            0.01975,
            0.04586,
            0.0483111,
        ],
    }
).set_index("terms")

comparison = (
    mulea_reference[["foreground", "eFDR"]]
    .rename(columns={"foreground": "mulea foreground", "eFDR": "mulea eFDR"})
    .join(
        efdr[["foreground", "padj"]].rename(
            columns={"foreground": "acore foreground", "padj": "acore eFDR"}
        )
    )
    .sort_values("mulea eFDR")
)
comparison

# %% [markdown]
# The same ten transcription factors are significant with `acore`, with the same
# number of foreground genes. The eFDR values agree up to the randomness of the
# permutations (and the fact that `acore` enforces a monotone eFDR).

# %%
set(efdr.index[efdr["rejected"]]) == set(mulea_reference.index)

# %%
ax = comparison.plot.scatter(
    x="mulea eFDR",
    y="acore eFDR",
    title="eFDR: mulea vs. acore",
    logx=False,
)
_ = ax.axline((0, 0), slope=1, color="grey", linestyle="--")

# %% [markdown]
# # Comparison with Benjamini-Hochberg
# The classic BH procedure on the same one-sided p-values gives a different number of
# significant transcription factors depending on which p-values are corrected: `acore`
# only tests annotations with at least `min_detected_in_set` foreground features (51
# here), whereas correcting over all 154 transcription factors of the ontology, where
# those without any gene in the target set have a p-value of 1, is more conservative.
# The eFDR on the other hand is based on the permuted data and not on the number of
# tests.
#
# If we would assign terms that do not have any foreground feature or to few detected
# in the set, then we would see the behaviour described by mulea w.r.t. the number
# of significant terms:
#
# ![eFDR normally more than BH](https://github.com/ELTEbioinformatics/mulea/raw/master/vignettes/Venn.png)

# %% tags=["hide-input"]
bh = acore.enrichment_analysis.run_regulation_enrichment(
    regulation_data=regulation_data,
    annotation=annotation.copy(),
    min_detected_in_set=1,
    correction="fdr_bh",
).set_index("terms")
n_terms_ontology = annotation["annotation"].nunique()
# terms without any gene in the foreground have a p-value of 1
pvalues_all_terms = bh["pvalue"].reindex(
    annotation["annotation"].unique(), fill_value=1.0
)
rejected_all_terms, _ = apply_pvalue_correction(
    pvalues_all_terms.to_numpy(), alpha=0.05, method="fdr_bh"
)
pd.Series(
    {
        "BH, tested terms only (acore default)": int(bh["rejected"].astype(bool).sum()),
        f"BH, all {n_terms_ontology} terms of the ontology (as in mulea)": int(
            rejected_all_terms.sum()
        ),
        "eFDR": int(efdr["rejected"].sum()),
    },
    name="significant at 0.05",
).to_frame()

# %% [markdown]
# # The eFDR step by step
# Here we compute the eFDR "by hand" for the same data to see what happens inside
# `correction="efdr"`. The idea in one sentence: *if the target set was just a random
# selection of genes, how many transcription factors would look as significant as the
# one we are looking at, compared to how many we observe?*
#
# Reference: [Turek et al. 2024, Sec. 2 on Implementation](https://pmc.ncbi.nlm.nih.gov/articles/PMC11490090/#Sec2)
#
# ## Step 0: which genes belong to which transcription factor
# A table of genes (rows) and transcription factors (columns) with `True` if the gene
# is regulated by the transcription factor. With it, the number of target genes
# regulated by each transcription factor is a column sum.

# %%
import numpy as np
from scipy import stats

from acore.multiple_testing import compute_efdr

genes = pd.Index(background_set)
tf_names = pd.Index(annotation["annotation"].unique())
# genes of the ontology which were not measured are not part of the test
measured = annotation[annotation["identifier"].isin(genes)]
membership = pd.DataFrame(False, index=genes, columns=tf_names)
membership.values[
    genes.get_indexer(measured["identifier"]),
    tf_names.get_indexer(measured["annotation"]),
] = True
membership.index.name = "gene"
membership.columns.name = "transcription factor"
n_genes, n_target = len(genes), len(target_set)
genes_per_tf = membership.sum()  # size of each regulon in the background
membership.iloc[4544:4550, :6]  # first 6 genes and transcription factors

# %% [markdown]
# ## Step 1: the observed p-values
# For each transcription factor we count how many of the 241 target genes it regulates
# and ask how surprising this is, if 241 genes were picked at random from all 7381
# genes. This is the one-sided Fisher's exact (hypergeometric) test. We only test
# transcription factors with at least one target gene (`min_detected_in_set=1`).


# %%
def enrichment_pvalues(n_hits: np.ndarray) -> np.ndarray:
    """P-value of seeing at least `n_hits` genes of a regulon in a random selection."""
    return stats.hypergeom.sf(n_hits - 1, n_genes, genes_per_tf.to_numpy(), n_target)


is_target = genes.isin(target_set)
observed = pd.DataFrame(
    {
        "regulon size": genes_per_tf,
        "target genes in regulon": membership[is_target].sum(),
    }
)
observed["expected by chance"] = observed["regulon size"] * n_target / n_genes
observed["pvalue"] = enrichment_pvalues(observed["target genes in regulon"].to_numpy())
observed = observed[observed["target genes in regulon"] >= 1].sort_values("pvalue")
print(f"{len(observed)} transcription factors with at least one target gene")
observed.head(6).round(5)

# %% [markdown]
# These are the same p-values as in the result of `acore` above:

# %%
np.allclose(observed["pvalue"].sort_index(), efdr["pvalue"].sort_index())

# %% [markdown]
# ## Step 2: one random selection of genes
# Now we pretend that the target set is a random selection of 241 genes. We repeat
# exactly the same tests for **all** transcription factors of the ontology, also those
# without hits (their p-value is 1). Of course, few transcription factors will look
# very significant by chance.

# %%
rng = np.random.default_rng(1)


def random_selection_pvalues() -> np.ndarray:
    """P-values of all transcription factors for a random selection of genes."""
    random_genes = rng.choice(n_genes, size=n_target, replace=False)
    n_hits = membership.to_numpy()[random_genes].sum(axis=0)
    return enrichment_pvalues(n_hits)


one_permutation = pd.Series(random_selection_pvalues(), index=tf_names)
one_permutation.sort_values().head(6).round(4)

# %% [markdown]
# We pick the p-value of one observed transcription factor, here the 3rd most
# significant one (`SoxS`), as a threshold $t$. In this single random selection
# the number of transcription factors with a p-value of at most $t$ is:

# %%
threshold = observed["pvalue"].iloc[2]
print(f"{observed.index[2]}: t = {threshold:.6f}")
(one_permutation <= threshold).sum()

# %% [markdown]
# ## Step 3: repeat many times
# One random selection is not enough. We repeat it 10000 times and get a table with a
# p-value for each transcription factor (columns) and random selection (rows).

# %%
n_permutations = 10_000
permuted = np.array([random_selection_pvalues() for _ in range(n_permutations)])
pd.DataFrame(permuted, columns=tf_names).head(6).round(4).rename_axis("permutation")

# %% [markdown]
# ## Step 4: how many false positives do we expect at a threshold?
# For the threshold $t$ of `SoxS`: count in each random selection how many
# transcription factors have a p-value of at most $t$. The distribution below shows
# that in most random selections there is none. The average over all random
# selections is the *expected number of false positives* at this threshold.

# %%
hits_at_threshold = (permuted <= threshold).sum(axis=1)
ax = (
    pd.Series(hits_at_threshold)
    .value_counts()
    .sort_index()
    .plot.bar(
        xlabel=f"transcription factors with p <= {threshold:.6f} in a random selection",
        ylabel=f"random selections (of {n_permutations})",
    )
)
expected_false_positives = hits_at_threshold.mean()
print(f"expected number of false positives at t: {expected_false_positives:.4f}")

# %% [markdown]
# ## Step 5: how many do we observe at the threshold?
# In the real data, `SoxS` is the third transcription factor with a p-value of at most
# $t$ (the number of observed transcription factors with $p \leq t$ is called
# $R_{obs}$ and the expected number in random selections $R_{exp}$). The estimated
# false discovery rate is the expected fraction of false positives among all calls
# at this threshold:
#
# $$eFDR(t) = \frac{R_{exp}(t)}{R_{obs}(t)}$$

# %%
n_observed = (observed["pvalue"] <= threshold).sum()
print(f"observed: {n_observed}, expected by chance: {expected_false_positives:.4f}")
print(f"eFDR: {expected_false_positives / n_observed:.5f}")

# %% [markdown]
# ## Step 6: do this for each transcription factor
# Use the p-value of each transcription factor as the threshold.

# %%
steps = pd.DataFrame({"pvalue": observed["pvalue"]})
steps["Rank_obs"] = [(observed["pvalue"] <= t).sum() for t in steps["pvalue"]]
steps["Rank_exp"] = [(permuted <= t).sum() / n_permutations for t in steps["pvalue"]]
steps["eFDR (raw)"] = steps["Rank_exp"] / steps["Rank_obs"]
steps.head(12).round(5)

# %% [markdown]
# ## Step 7: make it monotone and cap at 1
# A transcription factor with a smaller p-value should never get a *larger* eFDR than
# one with a larger p-value (see `MarA` and `IHF` in the table above).
# Therefore each eFDR is replaced by the smallest eFDR of all transcription factors with
# an equal or larger p-value (cumulative minimum from the largest p-value), and eFDRs
# are capped at 1. A transcription factor is significant if its eFDR is at most 0.05.

# %%
steps["eFDR"] = np.minimum(steps["eFDR (raw)"][::-1].cummin()[::-1], 1.0)
steps["significant"] = steps["eFDR"] <= 0.05
steps.head(12).round(5)

# %% [markdown]
# ## Check: the vectorized implementation gives the same result
# `acore.multiple_testing.compute_efdr` does the steps 4 to 7 for all transcription
# factors at once. Given the same p-values from the random selections, the result is
# identical to our step-by-step computation.

# %%
rejected, efdr_vectorized = compute_efdr(
    observed["pvalue"].to_numpy(), list(permuted), alpha=0.05
)
assert np.allclose(efdr_vectorized, steps["eFDR"])
assert (rejected == steps["significant"]).all()
print(f"identical eFDR for all {len(steps)} transcription factors")

# %% [markdown]
# The complete function (`run_regulation_enrichment(..., correction="efdr")`) uses
# other random selections than our walkthrough, so the eFDR values are
# similar but not identical, and the same transcription factors are significant:

# %%
manual_vs_function = steps[["eFDR"]].join(efdr["padj"].rename("eFDR (acore function)"))
print(
    "same significant transcription factors:",
    set(steps.index[steps["significant"]]) == set(efdr.index[efdr["rejected"]]),
)
manual_vs_function.head(12).round(5)

# %% [markdown]
# # Differences to be aware of
# - `acore` applies the multiple testing correction only to the annotations with at
#   least `min_detected_in_set` foreground features, mulea to all annotations of the
#   ontology. For the BH correction this changes the number of tests.
# - Both use random permutations, so the eFDR varies slightly between runs unless a
#   seed is set (and it differs between implementations as random numbers differ).
# - `acore` enforces that the eFDR does not decrease with increasing p-values, mulea
#   does not.
#
# The mulea code of the README for reference (R):
#
# ```r
# ora_model <- ora(gmt = tf_ontology_filtered,
#                  element_names = target_set,
#                  background_element_names = background_set,
#                  p_value_adjustment_method = "eFDR",
#                  number_of_permutations = 10000,
#                  nthreads = 2,
#                  random_seed = 1)
# ora_results <- run_test(ora_model)
# ```
