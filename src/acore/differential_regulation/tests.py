"""All the tests for differential regulation. Functions used in the user facing
function starting with `run_`.

"""

import re

import numpy as np
import pandas as pd
import pingouin as pg
import scipy.stats
from statsmodels.formula.api import ols

from acore.multiple_testing import (
    apply_pvalue_correction,
    apply_pvalue_permutation_fdrcorrection,
    get_max_permutations,
)


# njab.stats.groups_comparision.py (partly renamed functions)
def calc_means_between_groups(
    df: pd.DataFrame,
    boolean_array: pd.Series,
    event_names: tuple[str, str] = ("1", "0"),
) -> pd.DataFrame:
    """Mean comparison between groups"""
    sub = df.loc[boolean_array].describe().iloc[:3]
    sub["event"] = event_names[0]
    sub = sub.set_index("event", append=True).swaplevel()
    ret = sub
    sub = df.loc[~boolean_array].describe().iloc[:3]
    sub["event"] = event_names[1]
    sub = sub.set_index("event", append=True).swaplevel()
    ret = pd.concat([ret, sub])
    ret.columns.name = "variable"
    ret.index.names = ("event", "stats")
    return ret.T


def calc_ttest(
    df: pd.DataFrame, boolean_array: pd.Series, variables: list[str]
) -> pd.DataFrame:
    """Calculate t-test for each variable in `variables` between two groups defined
    by boolean array."""
    ret = []
    for var in variables:
        _ = pg.ttest(df.loc[boolean_array, var], df.loc[~boolean_array, var])
        ret.append(_)
    ret = pd.concat(ret)
    ret = ret.set_index(variables)
    ret.columns.name = "ttest"
    ret.columns = pd.MultiIndex.from_product(
        [["ttest"], ret.columns], names=("test", "var")
    )
    return ret


# end njab.stats.groups_comparision.py


def calculate_ttest(
    df,
    condition1,
    condition2,
    paired=False,
    is_logged=True,
    non_par=False,
    tail="two-sided",
    correction="auto",
    r=0.707,
):
    """
    Calculates the t-test for the means of independent samples belonging to two different
    groups using scipy.stats.ttest_ind_.
    
    .. _scipy.stats.ttest_ind: \
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html.

    :param df: pandas dataframe with groups and subjects as rows and protein identifier
               as column.
    :param str condition1: identifier of first group.
    :param str condition2: ientifier of second group.
    :param bool is_logged: data is logged transformed
    :param bool non_par: if True, normality and variance equality assumptions are checked
                         and non-parametric test Mann Whitney U test if not passed
    :return: Tuple with t-statistics, two-tailed p-value, mean of first group,
             mean of second group and logfc.

    Example::

        result = calculate_ttest(df, 'group1', 'group2')
    """
    t = None
    pvalue = np.nan
    group1 = df[[condition1]].values
    group2 = df[[condition2]].values

    mean1 = group1.mean()
    std1 = group1.std()
    mean2 = group2.mean()
    std2 = group2.std()
    if is_logged:
        fc = mean1 - mean2
    else:
        fc = mean1 / mean2

    test = "t-Test"
    if not non_par:
        result = pg.ttest(
            df[condition1],
            df[condition2],
            paired=paired,
            alternative=tail,
            correction=correction,
            r=r,
        )
    else:
        test = "Mann Whitney"
        result = pg.mwu(group1, group2, alternative=tail)

    if "T" in result.columns:
        t = result["T"].values[0]
    elif "U-val" in result.columns:
        t = result["U-val"].values[0]
    if "p-val" in result.columns:
        pvalue = result["p-val"].values[0]

    return (t, pvalue, mean1, mean2, std1, std2, fc, test)


def calculate_thsd(df, column, group="group", alpha=0.05, is_logged=True):
    """
    Pairwise Tukey-HSD posthoc test using pingouin.pairwise_tukey_.

    .. _pingouin.pairwise_tukey: \
    https://pingouin-stats.org/build/html/generated/pingouin.pairwise_tukey.html

    :param df: pandas dataframe with group and protein identifier as columns
    :param str column: column containing the protein identifier
    :param str group: column label containing the between factor
    :param float alpha: significance level
    :return: Pandas dataframe.

    Example::

        result = calculate_thsd(df, column='HBG2~P69892', group='group', alpha=0.05)
    """
    posthoc = None
    posthoc = pg.pairwise_tukey(data=df, dv=column, between=group)
    posthoc.columns = [
        "group1",
        "group2",
        "mean(group1)",
        "mean(group2)",
        "log2FC",
        "std_error",
        "t-statistics",
        "posthoc pvalue",
        "effsize",
    ]
    posthoc["efftype"] = "hedges"
    posthoc = complement_posthoc(posthoc, identifier=column, is_logged=is_logged)

    return posthoc


def calculate_pairwise_ttest(
    df, column, subject="subject", group="group", correction="none", is_logged=True
):
    """
    Performs pairwise t-test using pingouin, as a posthoc test,
    and calculates fold-changes using pingouin.pairwise_ttests_.
    
    .. _pingouin.pairwise_ttests: \
    https://pingouin-stats.org/build/html/generated/pingouin.pairwise_ttests.html.

    :param df: pandas dataframe with subject and group as rows and protein identifier as column.
    :param str column: column label containing the dependant variable
    :param str subject: column label containing subject identifiers
    :param str group: column label containing the between factor
    :param str correction: method used for testing and adjustment of p-values.
    :return: Pandas dataframe with means, standard deviations, test-statistics, 
             degrees of freedom and effect size columns.

    Example::

        result = calculate_pairwise_ttest(df,
                                          'protein a',
                                          subject='subject',
                                          group='group',
                                          correction='none'
                )
    """

    posthoc_columns = [
        "Contrast",
        "group1",
        "group2",
        "mean(group1)",
        "std(group1)",
        "mean(group2)",
        "std(group2)",
        "posthoc Paired",
        "posthoc Parametric",
        "posthoc T-Statistics",
        "posthoc dof",
        "posthoc tail",
        "posthoc pvalue",
        "posthoc BF10",
        "posthoc effsize",
    ]
    valid_cols = [
        "group1",
        "group2",
        "mean(group1)",
        "std(group1)",
        "mean(group2)",
        "std(group2)",
        "posthoc Paired",
        "posthoc Parametric",
        "posthoc T-Statistics",
        "posthoc dof",
        "posthoc tail",
        "posthoc pvalue",
        "posthoc BF10",
        "posthoc effsize",
    ]
    posthoc = df.pairwise_tests(
        dv=column,
        between=group,
        subject=subject,
        effsize="hedges",
        return_desc=True,
        padjust=correction,
    )
    posthoc.columns = posthoc_columns
    posthoc = posthoc[valid_cols]
    posthoc = complement_posthoc(posthoc, column, is_logged)
    posthoc["efftype"] = "hedges"

    return posthoc


def complement_posthoc(posthoc, identifier, is_logged):
    """
    Calculates fold-changes after posthoc test.

    :param posthoc: pandas dataframe from posthoc test. Should have at least columns
                    'mean(group1)' and 'mean(group2)'.
    :param str identifier: feature identifier.
    :return: Pandas dataframe with additional columns 'identifier', 'log2FC' and 'FC'.
    """
    posthoc["identifier"] = identifier
    if is_logged:
        posthoc["log2FC"] = posthoc["mean(group1)"] - posthoc["mean(group2)"]
        posthoc["FC"] = posthoc["log2FC"].apply(lambda x: np.power(2, x))
    else:
        posthoc["FC"] = posthoc["mean(group1)"] / posthoc["mean(group2)"]

    return posthoc


def calculate_anova(df, column, group="group"):
    """
    Calculates one-way ANOVA using pingouin.

    :param df: pandas dataframe with group as rows and protein identifier as column
    :param str column: name of the column in df to run ANOVA on
    :param str group: column with group identifiers
    :return: Tuple with t-statistics and p-value.
    """
    aov_result = pg.anova(data=df, dv=column, between=group)
    sel_cols = ["ddof1", "ddof2", "F", "p-unc"]
    df1, df2, t, pvalue = aov_result[sel_cols].values.tolist()[0]

    return (column, df1, df2, t, pvalue)


def calculate_ancova(data, column, group="group", covariates=[]):
    """
    Calculates one-way ANCOVA using pingouin.

    :param df: pandas dataframe with group as rows and protein identifier as column
    :param str column: name of the column in df to run ANOVA on
    :param str group: column with group identifiers
    :param list covariates: list of covariates (columns in df)
    :return: Tuple with column, F-statistics and p-value.
    """
    ancova_result = pg.ancova(data=data, dv=column, between=group, covar=covariates)
    t, df, pvalue = (
        ancova_result.loc[ancova_result["Source"] == group, ["F", "DF", "p-unc"]]
        .values.tolist()
        .pop()
    )

    return (column, df, df, t, pvalue)


def calculate_repeated_measures_anova(df, column, subject="subject", within="group"):
    """
    One-way and two-way repeated measures ANOVA using pingouin stats.

    :param df: pandas dataframe with samples as rows and protein identifier as column.
               Data must be in long-format for two-way repeated measures.
    :param str column: column label containing the dependant variable
    :param str subject: column label containing subject identifiers
    :param str within: column label containing the within factor
    :return: Tuple with protein identifier, t-statistics and p-value.

    Example::

        result = calculate_repeated_measures_anova(df,
                                                  'protein a',
                                                  subject='subject',
                                                  within='group'
                )
    """
    df1 = np.nan
    df2 = np.nan
    t = np.nan
    pvalue = np.nan
    try:
        aov_result = pg.rm_anova(
            data=df,
            dv=column,
            within=within,
            subject=subject,
            detailed=True,
            correction=True,
        )
        t, pvalue = aov_result.loc[0, ["F", "p-unc"]].values.tolist()
        df1, df2 = aov_result["DF"]
    except Exception as e:
        print(
            f"Repeated measurements Anova for column: {column} could not be calculated."
            f" Error {e}"
        )

    return (column, df1, df2, t, pvalue)


def calculate_mixed_anova(
    df, column, subject="subject", within="group", between="group2"
):
    """
    One-way and two-way repeated measures ANOVA using pingouin stats.

    :param df: pandas dataframe with samples as rows and protein identifier as column.
               Data must be in long-format for two-way repeated measures.
    :param str column: column label containing the dependant variable
    :param str subject: column label containing subject identifiers
    :param str within: column label containing the within factor
    :param str within: column label containing the between factor
    :return: Tuple with protein identifier, t-statistics and p-value.

    Example::

        result = calculate_mixed_anova(df,
                                       'protein a',
                                       subject='subject',
                                       within='group',
                                       between='group2'
                )
    """
    try:
        aov_result = pg.mixed_anova(
            data=df,
            dv=column,
            within=within,
            between=between,
            subject=subject,
            correction=True,
        )
        aov_result["identifier"] = column
    except Exception as e:
        print(f"Mixed Anova for column: {column} could not be calculated. Error {e}")

    return aov_result[["identifier", "DF1", "DF2", "F", "p-unc", "Source"]]


def pairwise_ttest_with_covariates(df, column, group, covariates, is_logged):
    """Pairwise t-test with covariates using statsmodels."""
    formula = f"Q('{column}') ~ C(Q('{group}'))"
    for c in covariates:
        formula += f" + Q('{c}')"
    model = ols(formula, data=df).fit()
    pw = model.t_test_pairwise(f"C(Q('{group}'))").result_frame
    pw = pw.reset_index()
    groups = "|".join(
        [re.escape(str(s)) for s in df[group].unique().tolist()]
    )  # ! this caused issues:
    regex = rf"({groups})\-({groups})"
    pw["group1"] = pw["index"].apply(lambda x: re.search(regex, x).group(2))
    pw["group2"] = pw["index"].apply(lambda x: re.search(regex, x).group(1))

    means = df.groupby(group)[column].mean().to_dict()
    stds = df.groupby(group)[column].std().to_dict()
    pw["mean(group1)"] = [means[g] for g in pw["group1"].tolist()]
    pw["mean(group2)"] = [means[g] for g in pw["group2"].tolist()]
    pw["std(group1)"] = [stds[g] for g in pw["group1"].tolist()]
    pw["std(group2)"] = [stds[g] for g in pw["group2"].tolist()]
    pw = pw.drop(["pvalue-hs", "reject-hs"], axis=1)
    pw = pw.rename(columns={"t": "posthoc T-Statistics", "P>|t|": "posthoc pvalue"})

    pw = pw[
        [
            "group1",
            "group2",
            "mean(group1)",
            "std(group1)",
            "mean(group2)",
            "std(group2)",
            "posthoc T-Statistics",
            "posthoc pvalue",
            "coef",
            "std err",
            "Conf. Int. Low",
            "Conf. Int. Upp.",
        ]
    ]
    pw = complement_posthoc(pw, column, is_logged)

    return pw


def format_anova_table(
    df,
    aov_results,
    pairwise_results,
    pairwise_cols,
    group,
    permutations,
    alpha,
    correction,
):
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
            if max_perm < permutations:
                permutations = max_perm
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
        scores["correction"] = "FDR correction BH"
        scores["padj"] = padj
        corrected = True

    res = pd.DataFrame(pairwise_results, columns=pairwise_cols).set_index("identifier")
    if not res.empty:
        res = res.join(scores[["F-statistics", "pvalue", "padj"]].astype("float"))
        res["correction"] = scores["correction"]
    else:
        res = scores
        res["log2FC"] = np.nan

    res = res.reset_index()
    res["rejected"] = res["padj"] < alpha

    if "posthoc pvalue" in res.columns:
        res["-log10 pvalue"] = [-np.log10(x) for x in res["posthoc pvalue"].values]
    else:
        res["-log10 pvalue"] = [-np.log10(x) for x in res["pvalue"].values]

    return res


def calculate_pvalue_from_tstats(tstat, dfn):
    """
    Calculate two-tailed p-values from T- or F-statistics.

    tstat: T/F distribution
    dfn: degrees of freedrom *n* (values) per protein (keys),
         i.e. number of obervations - number of groups (dict)
    """
    pval = scipy.stats.t.sf(np.abs(tstat), dfn) * 2

    return pval


def eta_squared(aov):
    """
    Calculates the effect size using Eta-squared.

    :param aov: pandas dataframe with anova results from statsmodels.
    :return: Pandas dataframe with additional Eta-squared column.
    """
    aov["eta_sq"] = "NaN"
    aov["eta_sq"] = aov[:-1]["sum_sq"] / sum(aov["sum_sq"])
    return aov


def omega_squared(aov):
    """
    Calculates the effect size using Omega-squared.

    :param aov: pandas dataframe with anova results from statsmodels.
    :return: Pandas dataframe with additional Omega-squared column.
    """
    mse = aov["sum_sq"][-1] / aov["df"][-1]
    aov["omega_sq"] = "NaN"
    aov["omega_sq"] = (aov[:-1]["sum_sq"] - (aov[:-1]["df"] * mse)) / (
        sum(aov["sum_sq"]) + mse
    )
    return aov
