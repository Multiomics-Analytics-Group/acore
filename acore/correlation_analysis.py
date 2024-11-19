import itertools

import numpy as np
import pandas as pd
import pingouin as pg
from scipy import stats
from scipy.special import betainc

import acore.utils as utils
from acore.multiple_testing import apply_pvalue_correction


def corr_lower_triangle(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """Compute the correlation matrix, returning only unique values (lower triangle).
    Passes kwargs to pandas.DataFrame.corr method.
    """
    corr_df = df.corr(**kwargs)
    lower_triangle = pd.DataFrame(np.tril(np.ones(corr_df.shape), -1)).astype(bool)
    lower_triangle.index, lower_triangle.columns = corr_df.index, corr_df.columns
    return corr_df.where(lower_triangle)


def calculate_correlations(x, y, method="pearson"):
    """
    Calculates a Spearman (nonparametric) or a Pearson (parametric) correlation coefficient and p-value to test for non-correlation.

    :param numpy.ndarray x: array 1
    :param numpy.ndarray y: array 2
    :param str method: chooses which kind of correlation method to run
    :return: Tuple with two floats, correlation coefficient and two-tailed p-value.

    Example::

        result = calculate_correlations(x, y, method='pearson')
    """
    if method == "pearson":
        coefficient, pvalue = stats.pearsonr(x, y)
    elif method == "spearman":
        coefficient, pvalue = stats.spearmanr(x, y)

    return (coefficient, pvalue)


def run_correlation(
    df,
    alpha=0.05,
    subject="subject",
    group="group",
    method="pearson",
    correction="fdr_bh",
):
    """
    This function calculates pairwise correlations for columns in dataframe, and returns it in the shape of a edge list with 'weight' as correlation score, and the ajusted p-values.

    :param df: pandas dataframe with samples as rows and features as columns.
    :param str subject: name of column containing subject identifiers.
    :param str group: name of column containing group identifiers.
    :param str method: method to use for correlation calculation ('pearson', 'spearman').
    :param float alpha: error rate. Values velow alpha are considered significant.
    :param str correction: type of correction see apply_pvalue_correction for methods
    :return: Pandas dataframe with columns: 'node1', 'node2', 'weight', 'padj' and 'rejected'.

    Example::

        result = run_correlation(df, alpha=0.05, subject='subject', group='group', method='pearson', correction='fdr_bh')
    """
    correlation = pd.DataFrame()
    # ToDo
    # The Repeated measurements correlation calculation is too time consuming so it only runs if
    # the number of features is less than 200
    if utils.check_is_paired(df, subject, group):
        if len(df[subject].unique()) > 2:
            if len(df.columns) < 200:
                correlation = run_rm_correlation(
                    df, alpha=alpha, subject=subject, correction=correction
                )
    else:
        df = df.dropna(axis=1)._get_numeric_data()
        if not df.empty:
            r, p = run_efficient_correlation(df, method=method)
            rdf = pd.DataFrame(r, index=df.columns, columns=df.columns)
            pdf = pd.DataFrame(p, index=df.columns, columns=df.columns)
            correlation = utils.convertToEdgeList(rdf, ["node1", "node2", "weight"])
            pvalues = utils.convertToEdgeList(pdf, ["node1", "node2", "pvalue"])
            correlation = pd.merge(correlation, pvalues, on=["node1", "node2"])

            rejected, padj = apply_pvalue_correction(
                correlation["pvalue"].tolist(), alpha=alpha, method=correction
            )
            correlation["padj"] = padj
            correlation["rejected"] = rejected
            correlation = correlation[correlation.rejected]
            correlation["pvalue"] = correlation["pvalue"].apply(
                lambda x: str(round(x, 5))
            )
            correlation["padj"] = correlation["padj"].apply(lambda x: str(round(x, 5)))

    return correlation


def run_multi_correlation(
    df_dict,
    alpha=0.05,
    subject="subject",
    on=["subject", "biological_sample"],
    group="group",
    method="pearson",
    correction="fdr_bh",
):
    """
    This function merges all input dataframes and calculates pairwise correlations for all columns.

    :param dict df_dict: dictionary of pandas dataframes with samples as rows and features as columns.
    :param str subject: name of the column containing subject identifiers.
    :param str group: name of the column containing group identifiers.
    :param list on: column names to join dataframes on (must be found in all dataframes).
    :param str method: method to use for correlation calculation ('pearson', 'spearman').
    :param float alpha: error rate. Values velow alpha are considered significant.
    :param str correction: type of correction see apply_pvalue_correction for methods
    :return: Pandas dataframe with columns: 'node1', 'node2', 'weight', 'padj' and 'rejected'.

    Example::

        result = run_multi_correlation(df_dict, alpha=0.05, subject='subject', on=['subject', 'biological_sample'] , group='group', method='pearson', correction='fdr_bh')
    """
    multidf = pd.DataFrame()
    correlation = None
    for dtype in df_dict:
        if multidf.empty:
            if isinstance(df_dict[dtype], pd.DataFrame):
                multidf = df_dict[dtype]
        else:
            if isinstance(df_dict[dtype], pd.DataFrame):
                multidf = pd.merge(multidf, df_dict[dtype], how="inner", on=on)
    if not multidf.empty:
        correlation = run_correlation(
            multidf,
            alpha=alpha,
            subject=subject,
            group=group,
            method=method,
            correction=correction,
        )

    return correlation


def calculate_rm_correlation(df, x, y, subject):
    """
    Computes correlation and p-values between two columns a and b in df.

    :param df: pandas dataframe with subjects as rows and two features and columns.
    :param str x: feature a name.
    :param str y: feature b name.
    :param subject: column name containing the covariate variable.
    :return: Tuple with values for: feature a, feature b, correlation, p-value and degrees of freedom.

    Example::

        result = calculate_rm_correlation(df, x='feature a', y='feature b', subject='subject')
    """
    result = pg.rm_corr(data=df, x=x, y=y, subject=subject)

    return (
        x,
        y,
        result["r"].values[0],
        result["pval"].values[0],
        result["dof"].values[0],
    )


def run_rm_correlation(df, alpha=0.05, subject="subject", correction="fdr_bh"):
    """
    Computes pairwise repeated measurements correlations for all columns in dataframe, and returns results as an edge list with 'weight' as correlation score, p-values, degrees of freedom and ajusted p-values.

    :param df: pandas dataframe with samples as rows and features as columns.
    :param str subject: name of column containing subject identifiers.
    :param float alpha: error rate. Values velow alpha are considered significant.
    :param str correction: type of correction type see apply_pvalue_correction for methods
    :return: Pandas dataframe with columns: 'node1', 'node2', 'weight', 'pvalue', 'dof', 'padj' and 'rejected'.

    Example::

        result = run_rm_correlation(df, alpha=0.05, subject='subject', correction='fdr_bh')
    """
    rows = []
    if not df.empty:
        df = df.set_index(subject)._get_numeric_data().dropna(axis=1)
        df.columns = df.columns.astype(str)
        combinations = itertools.combinations(df.columns, 2)
        df = df.reset_index()
        for x, y in combinations:
            row = [x, y]
            subset = df[[x, y, subject]]
            row.extend(pg.rm_corr(subset, x, y, subject).values.tolist()[0])
            rows.append(row)

        correlation = pd.DataFrame(
            rows,
            columns=["node1", "node2", "weight", "dof", "pvalue", "CI95%", "power"],
        )
        rejected, padj = apply_pvalue_correction(
            correlation["pvalue"].tolist(), alpha=alpha, method=correction
        )
        correlation["padj"] = padj
        correlation["rejected"] = rejected
        correlation = correlation[correlation.rejected]
        correlation["padj"] = correlation["padj"].apply(lambda x: str(round(x, 5)))

    return correlation


def run_efficient_correlation(data, method="pearson"):
    """
    Calculates pairwise correlations and returns lower triangle of the matrix with correlation values and p-values.

    :param data: pandas dataframe with samples as index and features as columns (numeric data only).
    :param str method: method to use for correlation calculation ('pearson', 'spearman').
    :return: Two numpy arrays: correlation and p-values.

    Example::

        result = run_efficient_correlation(data, method='pearson')
    """
    matrix = data.values
    if method == "pearson":
        r = np.corrcoef(matrix, rowvar=False)
    elif method == "spearman":
        r, p = stats.spearmanr(matrix, axis=0)

    diagonal = np.triu_indices(r.shape[0], 1)
    rf = r[diagonal]
    df = matrix.shape[1] - 2
    ts = rf * rf * (df / (1 - rf * rf))
    pf = betainc(0.5 * df, 0.5, df / (df + ts))
    p = np.zeros(shape=r.shape)
    p[np.triu_indices(p.shape[0], 1)] = pf
    p[np.tril_indices(p.shape[0], -1)] = pf
    p[np.diag_indices(p.shape[0])] = np.ones(p.shape[0])

    r[diagonal] = np.nan
    p[diagonal] = np.nan

    return r, p
