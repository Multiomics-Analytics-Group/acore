import re

import numpy as np
import pandas as pd
import pingouin as pg
import scipy.stats
import statsmodels.api as sm
from statsmodels.formula.api import ols

import acore.utils
from acore.multiple_testing import (
    apply_pvalue_correction,
    apply_pvalue_permutation_fdrcorrection,
    correct_pairwise_ttest,
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
    df: pd.DataFrame, boolean_array: pd.Series, vars: list[str]
) -> pd.DataFrame:
    """Calculate t-test for each variable in `vars` between two groups defined
    by boolean array."""
    ret = []
    for var in vars:
        _ = pg.ttest(df.loc[boolean_array, var], df.loc[~boolean_array, var])
        ret.append(_)
    ret = pd.concat(ret)
    ret = ret.set_index(vars)
    ret.columns.name = "ttest"
    ret.columns = pd.MultiIndex.from_product(
        [["ttest"], ret.columns], names=("test", "var")
    )
    return ret


def run_diff_analysis(
    df: pd.DataFrame,
    boolean_array: pd.Series,
    event_names: tuple[str, str] = ("1", "0"),
    ttest_vars=("alternative", "p-val", "cohen-d"),
) -> pd.DataFrame:
    """Differential analysis procedure between two groups. Calculaes
    mean per group and t-test for each variable in `vars` between two groups."""
    ret = calc_means_between_groups(
        df, boolean_array=boolean_array, event_names=event_names
    )
    ttests = calc_ttest(df, boolean_array=boolean_array, vars=ret.index)
    ret = ret.join(ttests.loc[:, pd.IndexSlice[:, ttest_vars]])
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
    groups. For more information visit
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


def calculate_THSD(df, column, group="group", alpha=0.05, is_logged=True):
    """
    Pairwise Tukey-HSD posthoc test using pingouin stats.
    For more information visit https://pingouin-stats.org/generated/pingouin.pairwise_tukey.html

    :param df: pandas dataframe with group and protein identifier as columns
    :param str column: column containing the protein identifier
    :param str group: column label containing the between factor
    :param float alpha: significance level
    :return: Pandas dataframe.

    Example::

        result = calculate_THSD(df, column='HBG2~P69892', group='group', alpha=0.05)
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
    Performs pairwise t-test using pingouin, as a posthoc test, and calculates fold-changes. For more information visit https://pingouin-stats.org/generated/pingouin.pairwise_ttests.html.

    :param df: pandas dataframe with subject and group as rows and protein identifier as column.
    :param str column: column label containing the dependant variable
    :param str subject: column label containing subject identifiers
    :param str group: column label containing the between factor
    :param str correction: method used for testing and adjustment of p-values.
    :return: Pandas dataframe with means, standard deviations, test-statistics, degrees of freedom and effect size columns.

    Example::

        result = calculate_pairwise_ttest(df, 'protein a', subject='subject', group='group', correction='none')
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
    posthoc = df.pairwise_ttests(
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

    :param posthoc: pandas dataframe from posthoc test. Should have at least columns 'mean(group1)' and 'mean(group2)'.
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
    df1, df2, t, pvalue = aov_result[["ddof1", "ddof2", "F", "p-unc"]].values.tolist()[
        0
    ]

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

    :param df: pandas dataframe with samples as rows and protein identifier as column. Data must be in long-format for two-way repeated measures.
    :param str column: column label containing the dependant variable
    :param str subject: column label containing subject identifiers
    :param str within: column label containing the within factor
    :return: Tuple with protein identifier, t-statistics and p-value.

    Example::

        result = calculate_repeated_measures_anova(df, 'protein a', subject='subject', within='group')
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
            "Repeated measurements Anova for column: {} could not be calculated."
            " Error {}".format(column, e)
        )

    return (column, df1, df2, t, pvalue)


def calculate_mixed_anova(
    df, column, subject="subject", within="group", between="group2"
):
    """
    One-way and two-way repeated measures ANOVA using pingouin stats.

    :param df: pandas dataframe with samples as rows and protein identifier as column. Data must be in long-format for two-way repeated measures.
    :param str column: column label containing the dependant variable
    :param str subject: column label containing subject identifiers
    :param str within: column label containing the within factor
    :param str within: column label containing the between factor
    :return: Tuple with protein identifier, t-statistics and p-value.

    Example::

        result = calculate_mixed_anova(df, 'protein a', subject='subject', within='group', between='group2')
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
        print(
            "Mixed Anova for column: {} could not be calculated. Error {}".format(
                column, e
            )
        )

    return aov_result[["identifier", "DF1", "DF2", "F", "p-unc", "Source"]]


def run_anova(
    df,
    alpha=0.05,
    drop_cols=["sample", "subject"],
    subject="subject",
    group="group",
    permutations=0,
    correction="fdr_bh",
    is_logged=True,
    non_par=False,
):
    """
    Performs statistical test for each protein in a dataset.
    Checks what type of data is the input (paired, unpaired or repeated measurements) and performs posthoc tests for multiclass data.
    Multiple hypothesis correction uses permutation-based if permutations>0 and Benjamini/Hochberg if permutations=0.

    :param df: pandas dataframe with samples as rows and protein identifiers as columns (with additional columns 'group', 'sample' and 'subject').
    :param str subject: column with subject identifiers
    :param str group: column with group identifiers
    :param list drop_cols: column labels to be dropped from the dataframe
    :param float alpha: error rate for multiple hypothesis correction
    :param int permutations: number of permutations used to estimate false discovery rates.
    :param bool non_par: if True, normality and variance equality assumptions are checked and non-parametric test Mann Whitney U test if not passed
    :return: Pandas dataframe with columns 'identifier', 'group1', 'group2', 'mean(group1)', 'mean(group2)', 'Log2FC', 'std_error', 'tail', 't-statistics', 'posthoc pvalue', 'effsize', 'efftype', 'FC', 'rejected', 'F-statistics', 'p-value', 'correction', '-log10 p-value', and 'method'.

    Example::

        result = run_anova(df, alpha=0.05, drop_cols=["sample",'subject'], subject='subject', group='group', permutations=50)
    """
    res = pd.DataFrame()
    if subject is not None and acore.utils.check_is_paired(df, subject, group):
        groups = df[group].unique()
        drop_cols = [d for d in drop_cols if d != subject]
        if len(groups) == 2:
            res = run_ttest(
                df,
                groups[0],
                groups[1],
                alpha=alpha,
                drop_cols=drop_cols,
                subject=subject,
                group=group,
                paired=True,
                correction=correction,
                permutations=permutations,
                is_logged=is_logged,
                non_par=non_par,
            )
        elif len(groups) > 2:
            res = run_repeated_measurements_anova(
                df,
                alpha=alpha,
                drop_cols=drop_cols,
                subject=subject,
                within=group,
                permutations=0,
                is_logged=is_logged,
            )
    elif len(df[group].unique()) == 2:
        groups = df[group].unique()
        drop_cols = [d for d in drop_cols if d != subject]
        res = run_ttest(
            df,
            groups[0],
            groups[1],
            alpha=alpha,
            drop_cols=drop_cols,
            subject=subject,
            group=group,
            paired=False,
            correction=correction,
            permutations=permutations,
            is_logged=is_logged,
            non_par=non_par,
        )
    elif len(df[group].unique()) > 2:
        df = df.drop(drop_cols, axis=1)
        aov_results = []
        pairwise_results = []
        for col in df.columns.drop(group).tolist():
            aov = calculate_anova(df[[group, col]], column=col, group=group)
            aov_results.append(aov)
            pairwise_result = calculate_pairwise_ttest(
                df[[group, col]],
                column=col,
                subject=subject,
                group=group,
                is_logged=is_logged,
            )
            pairwise_cols = pairwise_result.columns
            pairwise_results.extend(pairwise_result.values.tolist())
        df = df.set_index([group])
        res = format_anova_table(
            df,
            aov_results,
            pairwise_results,
            pairwise_cols,
            group,
            permutations,
            alpha,
            correction,
        )
        res["Method"] = "One-way anova"
        res = correct_pairwise_ttest(res, alpha, correction)

    return res


def run_ancova(
    df,
    covariates,
    alpha=0.05,
    drop_cols=["sample", "subject"],
    subject="subject",
    group="group",
    permutations=0,
    correction="fdr_bh",
    is_logged=True,
    non_par=False,
):
    """
    Performs statistical test for each protein in a dataset.
    Checks what type of data is the input (paired, unpaired or repeated measurements) and performs posthoc tests for multiclass data.
    Multiple hypothesis correction uses permutation-based if permutations>0 and Benjamini/Hochberg if permutations=0.

    :param df: pandas dataframe with samples as rows and protein identifiers and covariates as columns (with additional columns 'group', 'sample' and 'subject').
    :param list covariates: list of covariates to include in the model (column in df)
    :param str subject: column with subject identifiers
    :param str group: column with group identifiers
    :param list drop_cols: column labels to be dropped from the dataframe
    :param float alpha: error rate for multiple hypothesis correction
    :param int permutations: number of permutations used to estimate false discovery rates.
    :param bool non_par: if True, normality and variance equality assumptions are checked and non-parametric test Mann Whitney U test if not passed
    :return: Pandas dataframe with columns 'identifier', 'group1', 'group2', 'mean(group1)', 'mean(group2)', 'Log2FC', 'std_error', 'tail', 't-statistics', 'posthoc pvalue', 'effsize', 'efftype', 'FC', 'rejected', 'F-statistics', 'p-value', 'correction', '-log10 p-value', and 'method'.

    Example::

        result = run_ancova(df, covariates=['age'], alpha=0.05, drop_cols=["sample",'subject'], subject='subject', group='group', permutations=50)
    """
    df = df.drop(drop_cols, axis=1)
    for cova in covariates:
        if df[cova].dtype != np.number:
            df[cova] = pd.Categorical(df[cova])
            df[cova] = df[cova].cat.codes

    pairwise_results = []
    ancova_result = []
    for col in df.columns.tolist():
        if col not in covariates and col != group:
            ancova = calculate_ancova(
                df[[group, col] + covariates], col, group=group, covariates=covariates
            )
            ancova_result.append(ancova)
            pairwise_result = pairwise_ttest_with_covariates(
                df, column=col, group=group, covariates=covariates, is_logged=is_logged
            )
            pairwise_cols = pairwise_result.columns
            pairwise_results.extend(pairwise_result.values.tolist())
    df = df.set_index([group])
    res = format_anova_table(
        df,
        ancova_result,
        pairwise_results,
        pairwise_cols,
        group,
        permutations,
        alpha,
        correction,
    )
    res["Method"] = "One-way ancova"
    res = correct_pairwise_ttest(res, alpha, correction)

    return res


def pairwise_ttest_with_covariates(df, column, group, covariates, is_logged):
    formula = "Q('%s') ~ C(Q('%s'))" % (column, group)
    for c in covariates:
        formula += " + Q('%s')" % (c)
    model = ols(formula, data=df).fit()
    pw = model.t_test_pairwise("C(Q('%s'))" % (group)).result_frame
    pw = pw.reset_index()
    groups = "|".join([re.escape(s) for s in df[group].unique().tolist()])
    regex = r"({})\-({})".format(groups, groups)
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


def run_repeated_measurements_anova(
    df,
    alpha=0.05,
    drop_cols=["sample"],
    subject="subject",
    within="group",
    permutations=50,
    correction="fdr_bh",
    is_logged=True,
):
    """
    Performs repeated measurements anova and pairwise posthoc tests for each protein in dataframe.

    :param df: pandas dataframe with samples as rows and protein identifiers as columns (with additional columns 'group', 'sample' and 'subject').
    :param str subject: column with subject identifiers
    :param str within: column with within factor identifiers
    :param list drop_cols: column labels to be dropped from the dataframe
    :param float alpha: error rate for multiple hypothesis correction
    :param int permutations: number of permutations used to estimate false discovery rates
    :return: Pandas dataframe

    Example::

        result = run_repeated_measurements_anova(df, alpha=0.05, drop_cols=['sample'], subject='subject', within='group', permutations=50)
    """
    df = df.drop(drop_cols, axis=1).dropna(axis=1)
    aov_results = []
    pairwise_results = []
    index = [within, subject]
    for col in df.columns.drop(index).tolist():
        cols = index + [col]
        aov = calculate_repeated_measures_anova(
            df[cols], column=col, subject=subject, within=within
        )
        aov_results.append(aov)
        pairwise_result = calculate_pairwise_ttest(
            df[[within, subject, col]],
            subject=subject,
            column=col,
            group=within,
            is_logged=is_logged,
        )
        pairwise_cols = pairwise_result.columns
        pairwise_results.extend(pairwise_result.values.tolist())

    df = df.set_index([subject, within])
    res = format_anova_table(
        df,
        aov_results,
        pairwise_results,
        pairwise_cols,
        within,
        permutations,
        alpha,
        correction,
    )
    res["Method"] = "Repeated measurements anova"
    res = correct_pairwise_ttest(res, alpha, correction=correction)

    return res


def run_mixed_anova(
    df,
    alpha=0.05,
    drop_cols=["sample"],
    subject="subject",
    within="group",
    between="group2",
    correction="fdr_bh",
):
    """
    In statistics, a mixed-design analysis of variance model, also known as a split-plot ANOVA, is used to test
    for differences between two or more independent groups whilst subjecting participants to repeated measures.
    Thus, in a mixed-design ANOVA model, one factor (a fixed effects factor) is a between-subjects variable and the other
    (a random effects factor) is a within-subjects variable. Thus, overall, the model is a type of mixed-effects model.
    [source:https://en.wikipedia.org/wiki/Mixed-design_analysis_of_variance]

    :param df: pandas dataframe with samples as rows and protein identifiers as columns (with additional columns 'group', 'sample' and 'subject').
    :param str subject: column with subject identifiers
    :param str within: column with within factor identifiers
    :param str between: column with between factor identifiers
    :param list drop_cols: column labels to be dropped from the dataframe
    :param float alpha: error rate for multiple hypothesis correction
    :param int permutations: number of permutations used to estimate false discovery rates
    :return: Pandas dataframe

    Example::

        result = run_mixed_anova(df, alpha=0.05, drop_cols=['sample'], subject='subject', within='group', between='group2', permutations=50)
    """
    df = df.drop(drop_cols, axis=1).dropna(axis=1)
    aov_results = []
    index = [within, subject, between]
    for col in df.columns.drop(index).tolist():
        cols = index + [col]
        aov = calculate_mixed_anova(
            df[cols], column=col, subject=subject, within=within, between=between
        )
        aov_results.append(aov)

    res = pd.concat(aov_results)
    res = res[res["Source"] == "Interaction"]
    res = res[["identifier", "DF1", "DF2", "F", "p-unc"]]
    res.columns = ["identifier", "dfk", "dfn", "F-statistics", "pvalue"]
    _, padj = apply_pvalue_correction(
        res["pvalue"].tolist(), alpha=alpha, method=correction
    )
    res["correction"] = "FDR correction BH"
    res["padj"] = padj
    res["rejected"] = res["padj"] < alpha
    res["testing"] = "Interaction"
    res["within"] = ",".join(df[within].unique().tolist())
    res["between"] = ",".join(df[between].unique().tolist())

    return res


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
    Performs p-value correction (permutation-based and FDR) and converts pandas dataframe into final format.

    :param df: pandas dataframe with samples as rows and protein identifiers as columns (with additional columns 'group', 'sample' and 'subject').
    :param list[tuple] aov_results: list of tuples with anova results (one tuple per feature).
    :param list[dataframes] pairwise_results: list of pandas dataframes with posthoc tests results
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
            scores["correction"] = "permutation FDR ({} perm)".format(permutations)
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


def run_ttest(
    df,
    condition1,
    condition2,
    alpha=0.05,
    drop_cols=["sample"],
    subject="subject",
    group="group",
    paired=False,
    correction="fdr_bh",
    permutations=0,
    is_logged=True,
    non_par=False,
):
    """
    Runs t-test (paired/unpaired) for each protein in dataset and performs permutation-based (if permutations>0) or Benjamini/Hochberg (if permutations=0) multiple hypothesis correction.

    :param df: pandas dataframe with samples as rows and protein identifiers as columns (with additional columns 'group', 'sample' and 'subject').
    :param str condition1: first of two conditions of the independent variable
    :param str condition2: second of two conditions of the independent variable
    :param str subject: column with subject identifiers
    :param str group: column with group identifiers (independent variable)
    :param list drop_cols: column labels to be dropped from the dataframe
    :param bool paired: paired or unpaired samples
    :param str correction: method of pvalue correction see apply_pvalue_correction for methods
    :param float alpha: error rate for multiple hypothesis correction
    :param int permutations: number of permutations used to estimate false discovery rates.
    :param bool is_logged: data is log-transformed
    :param bool non_par: if True, normality and variance equality assumptions are checked and non-parametric test Mann Whitney U test if not passed
    :return: Pandas dataframe with columns 'identifier', 'group1', 'group2', 'mean(group1)', 'mean(group2)', 'std(group1)', 'std(group2)', 'Log2FC', 'FC', 'rejected', 'T-statistics', 'p-value', 'correction', '-log10 p-value', and 'method'.

    Example::

        result = run_ttest(df, condition1='group1', condition2='group2', alpha = 0.05, drop_cols=['sample'], subject='subject', group='group', paired=False, correction='fdr_bh', permutations=50)
    """
    columns = [
        "T-statistics",
        "pvalue",
        "mean_group1",
        "mean_group2",
        "std(group1)",
        "std(group2)",
        "log2FC",
        "test",
    ]
    df = df.set_index(group)
    df = df.drop(drop_cols, axis=1)
    method = "Unpaired t-test"
    if non_par:
        method = "Unpaired t-Test and Mann-Whitney U test"

    if paired:
        df = df.reset_index().set_index([group, subject])
        method = "Paired t-test"
    else:
        if subject is not None:
            df = df.drop([subject], axis=1)

    scores = df.T.apply(
        func=calculate_ttest,
        axis=1,
        result_type="expand",
        args=(condition1, condition2, paired, is_logged, non_par),
    )
    scores.columns = columns
    scores = scores.dropna(how="all")

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
            scores["correction"] = "permutation FDR ({} perm)".format(permutations)
            corrected = True

    if not corrected:
        rejected, padj = apply_pvalue_correction(
            scores["pvalue"].tolist(), alpha=alpha, method=correction
        )
        scores["correction"] = "FDR correction BH"
        scores["padj"] = padj
        scores["rejected"] = rejected
        corrected = True

    scores["group1"] = condition1
    scores["group2"] = condition2
    if is_logged:
        scores["FC"] = scores["log2FC"].apply(lambda x: np.power(2, x))
    else:
        scores = scores.rename(columns={"log2FC": "FC"})

    scores["-log10 pvalue"] = [
        -np.log10(x) if x != 0 else -np.log10(alpha) for x in scores["pvalue"].values
    ]
    scores["Method"] = method
    scores.index.name = "identifier"
    scores = scores.reset_index()

    return scores


def calculate_pvalue_from_tstats(tstat, dfn):
    """
    Calculate two-tailed p-values from T- or F-statistics.

    tstat: T/F distribution
    dfn: degrees of freedrom *n* (values) per protein (keys), i.e. number of obervations - number of groups (dict)
    """
    pval = scipy.stats.t.sf(np.abs(tstat), dfn) * 2

    return pval


def run_two_way_anova(
    df, drop_cols=["sample"], subject="subject", group=["group", "secondary_group"]
):
    """
    Run a 2-way ANOVA when data['secondary_group'] is not empty

    :param df: processed pandas dataframe with samples as rows, and proteins and groups as columns.
    :param list drop_cols: column names to drop from dataframe
    :param str subject: column name containing subject identifiers.
    :param list group: column names corresponding to independent variable groups
    :return: Two dataframes, anova results and residuals.

    Example::

        result = run_two_way_anova(data, drop_cols=['sample'], subject='subject', group=['group', 'secondary_group'])
    """
    data = df.copy()
    factorA, factorB = group
    data = data.set_index([subject] + group)
    data = data.drop(drop_cols, axis=1)
    data.columns = data.columns.str.replace(r"-", "_")

    aov_result = []
    residuals = {}
    for col in data.columns:
        model = ols(
            "{} ~ C({})*C({})".format(col, factorA, factorB),
            data[col].reset_index().sort_values(group, ascending=[True, False]),
        ).fit()
        aov_table = sm.stats.anova_lm(model, typ=2)
        eta_squared(aov_table)
        omega_squared(aov_table)
        for i in aov_table.index:
            if i != "Residual":
                t, p, eta, omega = aov_table.loc[
                    i, ["F", "PR(>F)", "eta_sq", "omega_sq"]
                ]
                protein = col.replace("_", "-")
                aov_result.append((protein, i, t, p, eta, omega))
        residuals[col] = model.resid

    anova_df = pd.DataFrame(
        aov_result,
        columns=[
            "identifier",
            "source",
            "F-statistics",
            "pvalue",
            "eta_sq",
            "omega_sq",
        ],
    )
    anova_df = anova_df.set_index("identifier")
    anova_df = anova_df.dropna(how="all")

    return anova_df, residuals


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
