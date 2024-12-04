"""Enrichment Analysis Module. Contains different functions to perform enrichment 
analysis.

Most things in this module are covered in https://www.youtube.com/watch?v=2NC1QOXmc5o
by Lars Juhl Jensen.
"""

from __future__ import annotations

import os
import re
import uuid

import gseapy as gp
import numpy as np
import pandas as pd
from scipy import stats

from acore.multiple_testing import apply_pvalue_correction

TYPE_COLS_MSG = """
columns: 'terms', 'identifiers', 'foreground',
    'background', foreground_pop, background_pop, 'pvalue', 'padj' and 'rejected'.
"""


def run_fisher(
    group1: list[int],
    group2: list[int],
    alternative: str = "two-sided",
) -> tuple[float, float]:
    """Run fisher's exact test on two groups using `scipy.stats.fisher_exact`.

    Example::

        # annotated   not-annotated
        # group1      a               b
        # group2      c               d


        odds, pvalue = stats.fisher_exact(group1=[a, b],
                                          group2 =[c, d]
                        )
    """

    odds, pvalue = stats.fisher_exact([group1, group2], alternative)

    return (odds, pvalue)


def run_kolmogorov_smirnov(dist1, dist2, alternative="two-sided"):
    """
    Compute the Kolmogorov-Smirnov statistic on 2 samples.
    See https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html

    :param list dist1: sequence of 1-D ndarray (first distribution to compare)
        drawn from a continuous distribution
    :param list dist2: sequence of 1-D ndarray (second distribution to compare)
        drawn from a continuous distribution
    :param str alternative: defines the alternative hypothesis (default is ‘two-sided’):
        * **'two-sided'**
        * **'less'**
        * **'greater'**
    :return: statistic float and KS statistic pvalue float Two-tailed p-value.

    Example::

        result = run_kolmogorov_smirnov(dist1, dist2, alternative='two-sided')

    """

    result = stats.ks_2samp(dist1, dist2, alternative=alternative, mode="auto")

    return result


def run_site_regulation_enrichment(
    regulation_data,
    annotation,
    identifier="identifier",
    groups=("group1", "group2"),
    annotation_col="annotation",
    rejected_col="rejected",
    group_col="group",
    method="fisher",
    regex="(\\w+~.+)_\\w\\d+\\-\\w+",  # ! add example to docstring of what kind of string this matches
    correction="fdr_bh",
):
    r"""
    This function runs a simple enrichment analysis for significantly
    regulated protein sites in a dataset.

    :param regulation_data: pandas.DataFrame resulting from differential
        regulation analysis.
    :param annotation: pandas.DataFrame with annotations for features
        (columns: 'annotation', 'identifier' (feature identifiers), and 'source').
    :param str identifier: name of the column from annotation containing
        feature identifiers.
    :param list groups: column names from regulation_data containing
        group identifiers.
    :param str annotation_col: name of the column from annotation containing
        annotation terms.
    :param str rejected_col: name of the column from regulation_data containing
        boolean for rejected null hypothesis.
    :param str group_col: column name for new column in annotation dataframe
        determining if feature belongs to foreground or background.
    :param str method: method used to compute enrichment
        (only 'fisher' is supported currently).
    :param str regex: how to extract the annotated identifier from the site identifier
    :return: pandas.DataFrame with columns: 'terms', 'identifiers', 'foreground',
        'background', foreground_pop, background_pop, 'pvalue', 'padj' and 'rejected'.

    Example::

        result = run_site_regulation_enrichment(regulation_data,
            annotation,
            identifier='identifier',
            groups=['group1', 'group2'],
            annotation_col='annotation',
            rejected_col='rejected',
            group_col='group',
            method='fisher',
            match="(\\w+~.+)_\\w\\d+\\-\\w+"
        )
    """
    result = pd.DataFrame()
    if regulation_data is not None:
        if not regulation_data.empty:
            new_ids = []
            for ident in regulation_data[identifier].tolist():
                match = re.search(regex, ident)
                if match is not None:
                    new_ids.append(match.group(1))
                else:
                    new_ids.append(ident)
            regulation_data[identifier] = new_ids
            regulation_data = regulation_data.drop_duplicates()
            result = run_regulation_enrichment(
                regulation_data,
                annotation,
                identifier,
                groups,
                annotation_col,
                rejected_col,
                group_col,
                method,
                correction,
            )

    return result


def run_up_down_regulation_enrichment(
    regulation_data,
    annotation,
    identifier="identifier",
    groups=("group1", "group2"),
    annotation_col="annotation",
    # rejected_col="rejected",
    group_col="group",
    method="fisher",
    correction="fdr_bh",
    alpha=0.05,
    lfc_cutoff=1,
):
    """
    This function runs a simple enrichment analysis for significantly regulated proteins
    distinguishing between up- and down-regulated.

    :param regulation_data: pandas.DataFrame resulting from differential regulation
        analysis (CKG's regulation table).
    :param annotation: pandas.DataFrame with annotations for features
        (columns: 'annotation', 'identifier' (feature identifiers), and 'source').
    :param str identifier: name of the column from annotation containing feature identifiers.
    :param list groups: column names from regulation_data containing group identifiers.
    :param str annotation_col: name of the column from annotation containing annotation terms.
    :param str rejected_col: name of the column from regulation_data containing boolean for
        rejected null hypothesis.
    :param str group_col: column name for new column in annotation dataframe determining
        if feature belongs to foreground or background.
    :param str method: method used to compute enrichment
        (only 'fisher' is supported currently).
    :param str correction: method to be used for multiple-testing correction
    :param float alpha: adjusted p-value cutoff to define significance
    :param float lfc_cutoff: log fold-change cutoff to define practical significance
    :return: pandas.DataFrame with columns: 'terms', 'identifiers', 'foreground',
        'background', 'pvalue', 'padj' and 'rejected'.

    Example::

        result = run_up_down_regulation_enrichment(
            regulation_data,
            annotation,
            identifier='identifier',
            groups=['group1',
            'group2'],
            annotation_col='annotation',
            rejected_col='rejected',
            group_col='group',
            method='fisher',
            correction='fdr_bh',
            alpha=0.05,
            lfc_cutoff=1,
        )
    """
    enrichment_results = {}
    for g1, g2 in regulation_data.groupby(groups).groups:
        df = regulation_data.groupby(groups).get_group((g1, g2))
        if "posthoc padj" in df:
            df["up_pairwise_regulation"] = (df["posthoc padj"] <= alpha) & (
                df["log2FC"] >= lfc_cutoff
            )
            df["down_pairwise_regulation"] = (df["posthoc padj"] <= alpha) & (
                df["log2FC"] <= -lfc_cutoff
            )
        else:
            df["up_pairwise_regulation"] = (df["padj"] <= alpha) & (
                df["log2FC"] >= lfc_cutoff
            )
            df["down_pairwise_regulation"] = (df["padj"] <= alpha) & (
                df["log2FC"] <= -lfc_cutoff
            )

        enrichment = run_regulation_enrichment(
            df,
            annotation,
            identifier=identifier,
            annotation_col=annotation_col,
            rejected_col="up_pairwise_regulation",
            group_col=group_col,
            method=method,
            correction=correction,
        )
        enrichment["direction"] = "upregulated"
        enrichment_results[g1 + "~" + g2] = enrichment
        enrichment = run_regulation_enrichment(
            df,
            annotation,
            identifier=identifier,
            annotation_col=annotation_col,
            rejected_col="down_pairwise_regulation",
            group_col=group_col,
            method=method,
            correction=correction,
        )
        enrichment["direction"] = "downregulated"
        enrichment_results[g1 + "~" + g2] = enrichment_results[g1 + "~" + g2].append(
            enrichment
        )

    return enrichment_results


# ! to move
def _annotate_features(
    features: pd.Series,
    in_foreground: set[str] | list[str],
    in_background: set[str] | list[str],
) -> pd.Series:
    """
    Annotate features as foreground or background based on their presence in the
    foreground and background lists.

    :param features: pandas.Series with features and their annotations.
    :param in_foreground: list of features identifiers in the foreground.
    :type in_foreground: set or list-like
    :param in_background: list of features identifiers in the background.
    :type in_background: set or list-like
    :return: pandas.Series containing 'foreground' or 'background'.
             missing values are preserved.

    Example::

        result = _annotate_features(features, in_foreground, in_background)
    """
    in_either_or = features.isin(in_foreground) | features.isin(in_background)
    res = (
        features.where(in_either_or, np.nan)
        .mask(features.isin(in_foreground), "foreground")
        .mask(features.isin(in_background), "background")
    )
    return res


def run_regulation_enrichment(
    regulation_data: pd.DataFrame,
    annotation: pd.DataFrame,
    identifier: str = "identifier",
    annotation_col: str = "annotation",
    rejected_col: str = "rejected",
    group_col: str = "group",
    method: str = "fisher",
    min_detected_in_set: int = 2,
    correction: str = "fdr_bh",
    correction_alpha: float = 0.05,
) -> pd.DataFrame:
    """
    This function runs a simple enrichment analysis for significantly regulated features
    in a dataset.

    :param regulation_data: pandas.DataFrame resulting from differential regulation analysis.
    :param annotation: pandas.DataFrame with annotations for features
        (columns: 'annotation', 'identifier' (feature identifiers), and 'source').
    :param str identifier: name of the column from annotation containing feature identifiers.
        It should also be present in `regulation_data`.
    :param str annotation_col: name of the column from annotation containing annotation terms.
    :param str rejected_col: name of the column from `regulation_data` containing boolean for
        rejected null hypothesis.
    :param str group_col: column name for new column in annotation dataframe determining
        if feature belongs to foreground or background.
    :param str method: method used to compute enrichment (only 'fisher' is supported currently).
    :param str correction: method to be used for multiple-testing correction
    :return: pandas.DataFrame with columns: 'terms', 'identifiers', 'foreground',
        'background', 'foreground_pop', 'background_pop', 'pvalue', 'padj' and 'rejected'.

    Example::

        result = run_regulation_enrichment(
            regulation_data,
            annotation,
            identifier='identifier',
            annotation_col='annotation',
            rejected_col='rejected',
            group_col='group',
            method='fisher',
            min_detected_in_set=2,
            correction='fdr_bh',
            correction_alpha=0.05,
         )
    """
    # ? can we remove NA features in that column?
    mask_rejected = regulation_data[rejected_col].astype(bool)
    foreground_list = regulation_data.loc[mask_rejected, identifier].unique()
    background_list = regulation_data.loc[~mask_rejected, identifier].unique()
    foreground_pop = len(foreground_list)
    background_pop = len(regulation_data[identifier].unique())
    # needs to allow for missing annotations
    annotation[group_col] = _annotate_features(
        features=annotation[identifier],
        in_foreground=foreground_list,
        in_background=background_list,
    )
    annotation = annotation.dropna(subset=[group_col])

    result = run_enrichment(
        annotation,
        foreground_id="foreground",
        background_id="background",
        foreground_pop=foreground_pop,
        background_pop=background_pop,
        annotation_col=annotation_col,
        group_col=group_col,
        identifier_col=identifier,
        method=method,
        correction=correction,
        min_detected_in_set=min_detected_in_set,
        correction_alpha=correction_alpha,
    )

    return result


def run_enrichment(
    data: pd.DataFrame,
    foreground_id: str,
    background_id: str,
    foreground_pop: int,
    background_pop: int,
    min_detected_in_set: int = 2,
    annotation_col: str = "annotation",
    group_col: str = "group",
    identifier_col: str = "identifier",
    method: str = "fisher",
    correction: str = "fdr_bh",
    correction_alpha: float = 0.05,
) -> pd.DataFrame:
    """
    Computes enrichment of the foreground relative to a given backgroung,
    using Fisher's exact test, and corrects for multiple hypothesis testing.

    :param data: pandas.DataFrame with annotations for dataset features
        (columns: 'annotation', 'identifier', 'group').
    :param str foreground_id: group identifier of features that belong to the foreground.
    :param str background_id: group identifier of features that belong to the background.
    :param int foreground_pop: number of features in the foreground.
    :param int background_pop: number of features in the background.
    :param str annotation_col: name of the column containing annotation terms.
    :param str group_col: name of column containing the group identifiers.
    :param str identifier_col: name of column containing dependent variables identifiers.
    :param str method: method used to compute enrichment (only 'fisher' is supported currently).
    :param str correction: method to be used for multiple-testing correction.
    :param float correction_alpha: adjusted p-value cutoff to define significance.
    :return: pandas.DataFrame with columns: annotation terms, features,
        number of foregroung/background features in each term,
        p-values and corrected p-values
        (columns: 'terms', 'identifiers', 'foreground',
        'background', 'pvalue', 'padj' and 'rejected').

    Example::

        result = run_enrichment(
            data,
            foreground='foreground',
            background='background',
            foreground_pop=len(foreground_list),
            background_pop=len(background_list),
            annotation_col='annotation',
            group_col='group',
            identifier_col='identifier',
            method='fisher',
         )
    """
    if method != "fisher":
        raise ValueError("Only Fisher's exact test is supported at the moment.")

    result = pd.DataFrame()
    terms = []
    ids = []
    pvalues = []
    fnum = []
    bnum = []
    countsdf = (
        data.groupby([annotation_col, group_col])
        .agg(["count"])[(identifier_col, "count")]
        .reset_index()
    )
    countsdf.columns = [annotation_col, group_col, "count"]
    for annotation in countsdf.loc[
        countsdf[group_col] == foreground_id, annotation_col
    ].unique():
        counts = countsdf[countsdf[annotation_col] == annotation]
        num_foreground = counts.loc[counts[group_col] == foreground_id, "count"].values
        num_background = counts.loc[counts[group_col] == background_id, "count"].values
        # ! counts should always be of length one count? squeeze?
        if len(num_foreground) == 1:
            num_foreground = num_foreground[0]
        if len(num_background) == 1:
            num_background = num_background[0]
        else:
            num_background = 0
        if num_foreground >= min_detected_in_set:
            _, pvalue = run_fisher(
                [num_foreground, foreground_pop - num_foreground],
                [num_background, background_pop - foreground_pop - num_background],
            )
            fnum.append(num_foreground)
            bnum.append(num_background)
            terms.append(annotation)
            pvalues.append(pvalue)
            ids.append(
                ",".join(
                    data.loc[
                        (data[annotation_col] == annotation)
                        & (data[group_col] == foreground_id),
                        identifier_col,
                    ]
                )
            )
    if len(pvalues) > 1:
        rejected, padj = apply_pvalue_correction(
            pvalues,
            alpha=correction_alpha,
            method=correction,
        )
        result = pd.DataFrame(
            {
                "terms": terms,
                "identifiers": ids,
                "foreground": fnum,
                "background": bnum,
                "foreground_pop": foreground_pop,
                "background_pop": background_pop,
                "pvalue": pvalues,
                "padj": padj,
                "rejected": rejected.astype(bool),
            }
        )
        result = result.sort_values(by="padj", ascending=True)

    return result


def run_ssgsea(
    data: pd.DataFrame,
    annotation: str,
    set_index: list[str],
    annotation_col: str = "an notation",
    identifier_col: str = "identifier",
    outdir: str = "tmp",
    min_size: int = 15,
    max_size: int = 500,
    scale: bool = False,
    permutations: int = 0,
) -> dict[str, pd.DataFrame]:
    """
    Project each sample within a data set onto a space of gene set enrichment scores using
    the single sample gene set enrichment analysis (ssGSEA) projection methodology
    described in Barbie et al., 2009:
    https://www.nature.com/articles/nature08460#Sec3 (search "Single Sample" GSEA).

    :param data: pandas.DataFrame with the quantified features (i.e. subject x proteins)
    :param annotation: pandas.DataFrame with the annotation to be used in the enrichment
        (i.e. CKG pathway annotation file)
    :param str annotation_col: name of the column containing annotation terms.
    :param str identifier_col: name of column containing dependent variables identifiers.
    :param list set_index: column/s to be used as index. Enrichment will be calculated
        for these values (i.e ["subject"] will return subjects x pathways matrix of
        enrichment scores)
    :param str out_dir: directory path where results will be stored
        (default None, tmp folder is used)
    :param int min_size: minimum number of features (i.e. proteins) in enriched terms
        (i.e. pathways)
    :param int max_size: maximum number of features (i.e. proteins) in enriched terms
        (i.e. pathways)
    :param bool scale: whether or not to scale the data
    :param int permutations: number of permutations used in the ssgsea analysis
    :return: dictionary with two dataframes: es - enrichment scores,
        and nes - normalized enrichment scores.

    Example::

        stproject = "P0000008"
        p = project.Project(
            stproject,
            datasets={},
            knowledge=None,
            report={},
            configuration_files=None,
        )
        p.build_project(False)
        p.generate_report()

        proteomics_dataset = p.get_dataset("proteomics")
        annotations = proteomics_dataset.get_dataframe("pathway annotation")
        processed = proteomics_dataset.get_dataframe('processed')

        result = run_ssgsea(
            processed,
            annotations,
            annotation_col='annotation',
            identifier_col='identifier',
            set_index=['group',
            'sample',
            'subject'],
            outdir=None,
            min_size=10,
            scale=False,
            permutations=0
        )
    """
    result = {}
    df = data.copy()
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Comine columns to create a unique name for each set (?)
    name = []
    index = data[set_index]
    for i, row in data[set_index].iterrows():
        name.append(
            "_".join(row[set_index].tolist())
        )  # this assumes strings as identifiers

    df["Name"] = name
    index.index = name
    df = df.drop(set_index, axis=1).set_index("Name").transpose()

    if annotation_col in annotation and identifier_col in annotation:
        grouped_annotations = (
            annotation.groupby(annotation_col)[identifier_col].apply(list).reset_index()
        )
        fid = uuid.uuid4()
        file_path = os.path.join(outdir, str(fid) + ".gmt")
        with open(file_path, "w", encoding="utf8") as out:
            for i, row in grouped_annotations.iterrows():
                out.write(
                    row[annotation_col]
                    + "\t"
                    + "\t".join(list(filter(None, row[identifier_col])))
                    + "\n"
                )
        enrichment = gp.ssgsea(
            data=df,
            gene_sets=str(file_path),
            outdir=outdir,
            min_size=min_size,
            max_size=max_size,
            scale=scale,
            permutation_num=permutations,
            no_plot=True,
            processes=1,
            seed=10,
            format="png",
        )

        enrichment_es = pd.DataFrame(enrichment.resultsOnSamples).transpose()
        enrichment_es = enrichment_es.join(index)
        enrichment_nes = enrichment.res2d.transpose()
        enrichment_nes = enrichment_nes.join(index)

        result = {"es": enrichment_es, "nes": enrichment_nes}

    return result
