"""Differential regulation module."""

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols

import acore.utils
from acore.multiple_testing import (
    apply_pvalue_correction,
    apply_pvalue_permutation_fdrcorrection,
    correct_pairwise_ttest,
    get_max_permutations,
)

from .tests import (  # calculate_thsd, complement_posthoc,
    calc_means_between_groups,
    calc_ttest,
    calculate_ancova,
    calculate_anova,
    calculate_mixed_anova,
    calculate_pairwise_ttest,
    calculate_repeated_measures_anova,
    calculate_ttest,
    eta_squared,
    format_anova_table,
    omega_squared,
    pairwise_ttest_with_covariates,
)

__all__ = [
    "run_anova",
    "run_ancova",
    "run_diff_analysis",
    "run_mixed_anova",
    "run_repeated_measurements_anova",
    "run_ttest",
    "run_two_way_anova",
]

# njab.stats.groups_comparision.py


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
    ttests = calc_ttest(df, boolean_array=boolean_array, variables=ret.index)
    ret = ret.join(ttests.loc[:, pd.IndexSlice[:, ttest_vars]])
    return ret


# ckg based:
def run_anova(
    df: pd.DataFrame,
    alpha: float = 0.05,
    drop_cols: list[str] = ["sample", "subject"],
    subject: str = "subject",
    group: str = "group",
    permutations: int = 0,
    correction: str = "fdr_bh",
    is_logged: bool = True,
    non_par: bool = False,
) -> pd.DataFrame:
    """
    Performs statistical test for each protein in a dataset.
    Checks what type of data is the input (paired, unpaired or repeated measurements) and
    performs posthoc tests for multiclass data.
    Multiple hypothesis correction uses permutation-based
    if permutations>0 and Benjamini/Hochberg if permutations=0.

    :param pd.DataFrame df: pandas dataframe with samples as rows and protein identifiers as columns
               (with additional columns 'group', 'sample' and 'subject').
    :param float alpha: error rate for multiple hypothesis correction
    :param list drop_cols: column labels to be dropped from the dataframe
    :param str subject: column with subject identifiers
    :param str group: column with group identifiers
    :param int permutations: number of permutations used to estimate false discovery rates.
    :param str correction: method of pvalue correction see apply_pvalue_correction for methods,
                            use methods available in acore.multiple_testing
    :param bool is_logged: whether data is log-transformed
    :param bool non_par: if True, normality and variance equality assumptions are checked
                         and non-parametric test Mann Whitney U test if not passed
    :return: Pandas dataframe with columns 'identifier', 'group1', 'group2',
        'mean(group1)', 'mean(group2)', 'Log2FC', 'std_error', 'tail', 't-statistics',
        'posthoc pvalue', 'effsize', 'efftype', 'FC', 'rejected', 'F-statistics', 'p-value',
        'correction', '-log10 p-value', and 'method'.

    Example::

        result = run_anova(df,
                           alpha=0.05,
                           drop_cols=["sample",'subject'],
                           subject='subject',
                           group='group',
                           permutations=50
                )
    """
    res = pd.DataFrame()
    if subject is not None and acore.utils.check_is_paired(df, subject, group):
        paired = True
    else:
        paired = False

    if len(df[group].unique()) == 2:
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
            paired=paired,
            correction=correction,
            permutations=permutations,
            is_logged=is_logged,
            non_par=non_par,
        )
    elif len(df[group].unique()) > 2:
        if paired:
            res = run_repeated_measurements_anova(
                df,
                alpha=alpha,
                drop_cols=drop_cols,
                subject=subject,
                within=group,
                permutations=0,
                is_logged=is_logged,
            )
        else:
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
    else:
        raise ValueError("Number of groups must be greater than 1")

    return res


def run_ancova(
    df: pd.DataFrame,
    covariates: list[str],
    alpha: float = 0.05,
    drop_cols: list[str] = ["sample", "subject"],
    subject: str = "subject",
    group: str = "group",
    permutations: int = 0,
    correction: str = "fdr_bh",
    is_logged: bool = True,
    non_par: bool = False,
) -> pd.DataFrame:
    """
    Performs statistical test for each protein in a dataset.
    Checks what type of data is the input (paired, unpaired or repeated measurements)
    and performs posthoc tests for multiclass data.
    Multiple hypothesis correction uses permutation-based
    if permutations>0 and Benjamini/Hochberg if permutations=0.

    :param pd.DataFrame df: Pandas DataFrame with samples as rows and protein identifiers and
               covariates as columns (with additional columns 'group', 'sample' and 'subject').
    :param list covariates: list of covariates to include in the model (column in df)
    :param float alpha: error rate for multiple hypothesis correction
    :param list drop_cols: column labels to be dropped from the DataFrame
    :param str subject: column with subject identifiers
    :param str group: column with group identifiers
    :param int permutations: number of permutations used to estimate false discovery rates.
    :param str correction: method of pvalue correction see apply_pvalue_correction for methods,
                           use methods available in acore.multiple_testing
    :param bool is_logged: whether data is log-transformed
    :param bool non_par: if True, normality and variance equality assumptions are checked
                         and non-parametric test Mann Whitney U test if not passed
    :return: Pandas DataFrame with columns 'identifier', 'group1', 'group2',
        'mean(group1)', 'mean(group2)', 'Log2FC', 'std_error', 'tail', 't-statistics',
        'posthoc pvalue', 'effsize', 'efftype', 'FC', 'rejected', 'F-statistics', 'p-value',
        'correction', '-log10 p-value', and 'method'.

    Example::

        result = run_ancova(df,
                            covariates=['age'],
                            alpha=0.05,
                            drop_cols=["sample",'subject'],
                            subject='subject',
                            group='group',
                            permutations=50
                )
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


def run_repeated_measurements_anova(
    df,
    alpha=0.05,
    drop_cols=["sample"],
    subject="subject",
    within="group",
    permutations=50,
    correction="fdr_bh",
    is_logged=True,
) -> pd.DataFrame:
    """
    Performs repeated measurements anova and pairwise posthoc tests for each protein in dataframe.

    :param pd.DataFrame df: Pandas DataFrame with samples as rows and protein identifiers as columns
               (with additional columns 'group', 'sample' and 'subject').
    :param float alpha: error rate for multiple hypothesis correction
    :param list drop_cols: column labels to be dropped from the DataFrame
    :param str subject: column with subject identifiers
    :param str within: column with within factor identifiers
    :param int permutations: number of permutations used to estimate false discovery rates
    :param str correction: method of pvalue correction see apply_pvalue_correction for methods,
                           use methods available in acore.multiple_testing
    :param bool is_logged: whether data is log-transformed
    :return: Pandas DataFrame

    Example::

        result = run_repeated_measurements_anova(df,
                                                 alpha=0.05,
                                                 drop_cols=['sample'],
                                                 subject='subject',
                                                 within='group',
                                                 permutations=50
                )
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
    In statistics, a mixed-design analysis of variance model, also known as a split-plot
    ANOVA, is used to test
    for differences between two or more independent groups whilst subjecting participants
    to repeated measures.
    Thus, in a mixed-design ANOVA model, one factor (a fixed effects factor) is a
    between-subjects variable and the other
    (a random effects factor) is a within-subjects variable. Thus, overall, the model is a
    type of mixed-effects model (source_)

    .. _source: https://en.wikipedia.org/wiki/Mixed-design_analysis_of_variance

    :param pd.DataFrame df: Pandas DataFrame with samples as rows and protein identifiers as columns
               (with additional columns 'group', 'sample' and 'subject').
    :param float alpha: error rate for multiple hypothesis correction
    :param list drop_cols: column labels to be dropped from the DataFrame
    :param str subject: column with subject identifiers
    :param str within: column with within factor identifiers
    :param str between: column with between factor identifiers
    :param str correction: method of pvalue correction see apply_pvalue_correction for methods,
                           use methods available in acore.multiple_testing
    :return: Pandas DataFrame
    :rtype: pd.DataFrame

    Example::

        result = run_mixed_anova(df,
                                 alpha=0.05,
                                 drop_cols=['sample'],
                                 subject='subject',
                                 within='group',
                                 between='group2',
                )
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
    Runs t-test (paired/unpaired) for each protein in dataset and performs
    permutation-based (if permutations>0) or Benjamini/Hochberg (if permutations=0)
    multiple hypothesis correction.

    :param pd.DataFrame df: Pandas DataFrame with samples as rows and protein identifiers as columns
               (with additional columns 'group', 'sample' and 'subject').
    :param str condition1: first of two conditions of the independent variable
    :param str condition2: second of two conditions of the independent variable
    :param float alpha: error rate for multiple hypothesis correction
    :param list drop_cols: column labels to be dropped from the DataFrame
    :param str subject: column with subject identifiers
    :param str group: column with group identifiers (independent variable)
    :param bool paired: paired or unpaired samples
    :param str correction: method of pvalue correction see apply_pvalue_correction for methods
    :param int permutations: number of permutations used to estimate false discovery rates.
    :param bool is_logged: data is log-transformed
    :param bool non_par: if True, normality and variance equality assumptions are checked
                         and non-parametric test Mann Whitney U test if not passed
    :return: Pandas DataFrame with columns 'identifier', 'group1', 'group2',
        'mean(group1)', 'mean(group2)', 'std(group1)', 'std(group2)', 'Log2FC', 'FC',
        'rejected', 'T-statistics', 'p-value', 'correction', '-log10 p-value', and 'method'.

    Example::

        result = run_ttest(df,
                           condition1='group1',
                           condition2='group2',
                           alpha = 0.05,
                           drop_cols=['sample'],
                           subject='subject',
                           group='group',
                           paired=False,
                           correction='fdr_bh',
                           permutations=50
                )
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
            scores["correction"] = f"permutation FDR ({permutations} perm)"
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


def run_two_way_anova(
    df,
    drop_cols=["sample"],
    subject="subject",
    group=["group", "secondary_group"],
):
    """
    Run a 2-way ANOVA when data['secondary_group'] is not empty

    :param pd.DataFrame df: processed pandas DataFrame with samples as rows,
               and proteins and groups as columns.
    :param list drop_cols: column names to drop from DataFrame
    :param str subject: column name containing subject identifiers.
    :param list group: column names corresponding to independent variable groups
    :return: Two DataFrames, anova results and residuals.

    Example::

        result = run_two_way_anova(data,
                                   drop_cols=['sample'],
                                   subject='subject',
                                   group=['group', 'secondary_group']
                )
    """
    data = df.copy()
    factor_a, factor_b = group
    data = data.set_index([subject] + group)
    data = data.drop(drop_cols, axis=1)
    data.columns = data.columns.str.replace(r"-", "_")

    aov_result = []
    residuals = {}
    for col in data.columns:
        model = ols(
            f"{col} ~ C({factor_a})*C({factor_b})",
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
