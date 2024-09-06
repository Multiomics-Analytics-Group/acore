import numpy as np
import pandas as pd
import pingouin as pg
from scipy.special import factorial
from sklearn.utils import shuffle
from statsmodels.stats import multitest


def apply_pvalue_correction(pvalues, alpha=0.05, method="bonferroni"):
    """
    Performs p-value correction using the specified method. For more information visit https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html.

    :param numpy.ndarray pvalues: et of p-values of the individual tests.
    :param float alpha: error rate.
    :param str method: method of p-value correction:
        - bonferroni : one-step correction
        - sidak : one-step correction
        - holm-sidak : step down method using Sidak adjustments
        - holm : step-down method using Bonferroni adjustments
        - simes-hochberg : step-up method (independent)
        - hommel : closed method based on Simes tests (non-negative)
        - fdr_bh : Benjamini/Hochberg (non-negative)
        - fdr_by : Benjamini/Yekutieli (negative)
        - fdr_tsbh : two stage fdr correction (non-negative)
        - fdr_tsbky : two stage fdr correction (non-negative)
    :return: Tuple with two arrays, boolen for rejecting H0 hypothesis and float for adjusted p-value.

    Exmaple::

        result = apply_pvalue_correction(pvalues, alpha=0.05, method='bonferroni')
    """
    p = np.array(pvalues)
    mask = np.isfinite(p)
    pval_corrected = np.full(p.shape, np.nan)
    pval_corrected[mask] = multitest.multipletests(p[mask], alpha, method)[1]
    rejected = [p <= alpha for p in pval_corrected]

    return (rejected, pval_corrected.tolist())


def apply_pvalue_fdrcorrection(pvalues, alpha=0.05, method="indep"):
    """
    Performs p-value correction for false discovery rate. For more information visit https://www.statsmodels.org/devel/generated/statsmodels.stats.multitest.fdrcorrection.html.

    :param numpy.ndarray pvalues: et of p-values of the individual tests.
    :param float alpha: error rate.
    :param str method: method of p-value correction ('indep', 'negcorr').
    :return: Tuple with two arrays, boolen for rejecting H0 hypothesis and float for adjusted p-value.

    Exmaple::

        result = apply_pvalue_fdrcorrection(pvalues, alpha=0.05, method='indep')
    """
    rejected, padj = multitest.fdrcorrection(pvalues, alpha, method)

    return (rejected, padj)


def apply_pvalue_twostage_fdrcorrection(pvalues, alpha=0.05, method="bh"):
    """
    Iterated two stage linear step-up procedure with estimation of number of true hypotheses. For more information visit https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.fdrcorrection_twostage.html.

    :param numpy.ndarray pvalues: et of p-values of the individual tests.
    :param float alpha: error rate.
    :param str method: method of p-value correction ('bky', 'bh').
    :return: Tuple with two arrays, boolen for rejecting H0 hypothesis and float for adjusted p-value.

    Exmaple::

        result = apply_pvalue_twostage_fdrcorrection(pvalues, alpha=0.05, method='bh')
    """
    rejected, padj, num_hyp, alpha_stages = multitest.fdrcorrection_twostage(
        pvalues, alpha, method
    )

    return (rejected, padj)


def apply_pvalue_permutation_fdrcorrection(
    df, observed_pvalues, group, alpha=0.05, permutations=50
):
    """
    This function applies multiple hypothesis testing correction using a permutation-based false discovery rate approach.

    :param df: pandas dataframe with samples as rows and features as columns.
    :param oberved_pvalues: pandas Series with p-values calculated on the originally measured data.
    :param str group: name of the column containing group identifiers.
    :param float alpha: error rate. Values velow alpha are considered significant.
    :param int permutations: number of permutations to be applied.
    :return: Pandas dataframe with adjusted p-values and rejected columns.

    Example::

        result = apply_pvalue_permutation_fdrcorrection(df, observed_pvalues, group='group', alpha=0.05, permutations=50)
    """
    i = permutations
    df_index = df.index.values
    df_columns = df.columns.values
    seen = [str([df_index] + df_columns.tolist())]
    rand_pvalues = []
    while i > 0:
        df_index = shuffle(df_index)
        df_columns = shuffle(df_columns)
        df_random = df.reset_index(drop=True)
        df_random.index = df_index
        df_random.index.name = group
        df_random.columns = df_columns
        rand_index = str([df_random.index] + df_columns.tolist())
        if rand_index not in seen:
            seen.append(rand_index)
            df_random = df_random.reset_index()
            for col in df_random.columns.drop(group):
                rand_pvalues.append(
                    calculate_anova(df_random, column=col, group=group)[-1]
                )
            i -= 1
    rand_pvalues = np.array(rand_pvalues)

    qvalues = []
    for i, row in observed_pvalues.to_frame():
        qvalues.append(
            get_counts_permutation_fdr(
                row["pvalue"], rand_pvalues, df["pvalue"], permutations, alpha
            )
            + (i,)
        )

    qvalues = pd.DataFrame(
        qvalues, columns=["padj", "rejected", "identifier"]
    ).set_index("identifier")

    return qvalues


def calculate_anova(df, column, group="group"):
    """
    Calculates one-way ANOVA using pingouin.

    :param df: pandas dataframe with group as rows and protein identifier as column
    :param str column: name of the column in df to run ANOVA on
    :param str group: column with group identifiers
    :return: Tuple with t-statistics and p-value.
    """
    aov_result = pg.anova(data=df, dv=column, between=group)
    df1, df2, t, pvalue = aov_result[["ddof1", "ddof2", "F", "p-unc"]].values.tolist()[
        0
    ]

    return (column, df1, df2, t, pvalue)


def get_counts_permutation_fdr(value, random, observed, n, alpha):
    """
    Calculates local FDR values (q-values) by computing the fraction of accepted hits from the permuted data over accepted hits from the measured data normalized by the total number of permutations.

    :param float value: computed p-value on measured data for a feature.
    :param numpy.ndarray random: p-values computed on the permuted data.
    :param observed: pandas Series with p-values calculated on the originally measured data.
    :param int n: number of permutations to be applied.
    :param float alpha: error rate. Values velow alpha are considered significant.
    :return: Tuple with q-value and boolean for H0 rejected.

    Example::

        result = get_counts_permutation_fdr(value, random, observed, n=250, alpha=0.05)
    """
    a = random[random <= value].shape[0] + 0.0000000000001  # Offset in case of a = 0.0
    b = (observed <= value).sum()
    qvalue = a / b / float(n)

    return (qvalue, qvalue <= alpha)


def get_max_permutations(df, group="group"):
    """
    Get maximum number of permutations according to number of samples.

    :param df: pandas dataframe with samples as rows and protein identifiers as columns
    :param str group: column with group identifiers
    :return: Maximum number of permutations.
    :rtype: int
    """
    num_groups = len(list(df.index))
    num_per_group = df.groupby(group).size().tolist()
    max_perm = factorial(num_groups) / np.prod(factorial(np.array(num_per_group)))

    return max_perm


def correct_pairwise_ttest(df, alpha, correction="fdr_bh"):
    posthoc_df = pd.DataFrame()
    if "group1" in df and "group2" in df and "posthoc pvalue" in df:
        for comparison in df.groupby(["group1", "group2"]).groups:
            index = df.groupby(["group1", "group2"]).groups.get(comparison)
            posthoc_pvalues = df.loc[index, "posthoc pvalue"].tolist()
            _, posthoc_padj = apply_pvalue_correction(
                posthoc_pvalues, alpha=alpha, method=correction
            )
            if posthoc_df.empty:
                posthoc_df = pd.DataFrame(
                    {"index": index, "posthoc padj": posthoc_padj}
                )
            else:
                posthoc_df = posthoc_df.append(
                    pd.DataFrame({"index": index, "posthoc padj": posthoc_padj})
                )
        posthoc_df = posthoc_df.set_index("index")
        df = df.join(posthoc_df)

    return df
