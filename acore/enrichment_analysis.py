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
    """Run fisher's exact test on two groups using `scipy.stats.fisher_exact`_.

    .. _scipy.stats.fisher_exact: https://docs.scipy.org/doc/scipy/reference/generated/\
scipy.stats.fisher_exact.html

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
    See `scipy.stats.ks_2samp`_

    .. _scipy.stats.ks_2samp: https://docs.scipy.org/doc/scipy/reference/generated/\
scipy.stats.ks_2samp.html


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


# ! undocumented for now (find usage example)
def run_site_regulation_enrichment(
    regulation_data: pd.DataFrame,
    annotation: pd.DataFrame,
    identifier: str = "identifier",
    groups: list[str] = ("group1", "group2"),
    annotation_col: str = "annotation",
    rejected_col: str = "rejected",
    group_col: str = "group",
    method: str = "fisher",
    regex: str = "(\\w+~.+)_\\w\\d+\\-\\w+",
    correction: str = "fdr_bh",
    remove_duplicates: bool = False,
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

    :raises ValueError: if regulation_data is `None` or empty.

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
    if regulation_data is None or regulation_data.empty:
        raise ValueError("regulation_data is empty")

    new_ids = []
    # find any identifiers with a PTM and save only prot+gene identifer
    for ident in regulation_data[identifier].tolist():
        match = re.search(regex, ident)
        if match is not None:
            new_ids.append(
                match.group(1)
            )  # removes the PTM extension of the identifier of CKG
        else:
            new_ids.append(ident)
    # so this is normalizing the identifiers to ignore the PTM extension
    regulation_data[identifier] = new_ids  # matches are used as identifiers
    if remove_duplicates:
        regulation_data = regulation_data.drop_duplicates(subset=[identifier])
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
    regulation_data: pd.DataFrame,
    annotation: pd.DataFrame,
    identifier: str = "identifier",
    groups: list[str] = ("group1", "group2"),
    annotation_col: str = "annotation",
    # rejected_col:str="rejected",
    group_col: str = "group",
    method: str = "fisher",
    min_detected_in_set: int = 2,
    correction: str = "fdr_bh",
    correction_alpha: float = 0.05,
    lfc_cutoff: float = 1,
) -> pd.DataFrame:
    """
    This function runs a simple enrichment analysis for significantly regulated proteins
    distinguishing between up- and down-regulated.

    :param pandas.DataFrame regulation_data: pandas.DataFrame resulting from differential regulation
        analysis (CKG's regulation table).
    :param pandas.DataFrame annotation: pandas.DataFrame with annotations for features
        (columns: 'annotation', 'identifier' (feature identifiers), and 'source').
    :param str identifier: name of the column from annotation containing feature identifiers.
    :param list[str] groups: column names from regulation_data containing group identifiers.
            See `pandas.DataFrame.groupby`_ for more information.
            
            .. _pandas.DataFrame.groupby: https://pandas.pydata.org/pandas-docs/stable/\
reference/api/pandas.DataFrame.groupby.html
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
        'background', 'pvalue', 'padj', 'rejected', 'direction' and 'comparison'.

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
    if isinstance(groups, str):
        groups = [groups]
    if isinstance(groups, tuple):
        groups = list(groups)
    if len(groups) != 2:
        raise ValueError("groups should contains exactly two columns.")

    ret = list()
    # In case of multiple comparisons this is used to get all possible combinations
    for g1, g2 in regulation_data.groupby(groups).groups:

        df = regulation_data.groupby(groups).get_group((g1, g2))

        padj_name = "padj"
        if "posthoc padj" in df:
            padj_name = "posthoc padj"

        df["up_pairwise_regulation"] = (df[padj_name] <= correction_alpha) & (
            df["log2FC"] >= lfc_cutoff
        )
        df["down_pairwise_regulation"] = (df[padj_name] <= correction_alpha) & (
            df["log2FC"] <= -lfc_cutoff
        )
        comparison_tag = g1 + "~" + g2

        for rej_col, direction in zip(
            ("up_pairwise_regulation", "down_pairwise_regulation"),
            ("upregulated", "downregulated"),
        ):
            _enrichment = run_regulation_enrichment(
                df,
                annotation,
                identifier=identifier,
                annotation_col=annotation_col,
                rejected_col=rej_col,
                group_col=group_col,
                method=method,
                min_detected_in_set=min_detected_in_set,
                correction=correction,
                correction_alpha=correction_alpha,
            )
            _enrichment["direction"] = direction
            _enrichment["comparison"] = comparison_tag
            ret.append(_enrichment)

    ret = pd.concat(ret)

    return ret


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
    if len(pvalues) >= 1:
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
    set_index: list[str] = None,
    annotation_col: str = "annotation",
    identifier_col: str = "identifier",
    outdir: str = "tmp",
    min_size: int = 15,
    max_size: int = 500,
    scale: bool = False,
    permutations: int = 0,
) -> pd.DataFrame:
    """
    Project each sample within a data set onto a space of gene set enrichment scores using
    the single sample gene set enrichment analysis (ssGSEA) projection methodology
    described in Barbie et al., 2009:
    https://www.nature.com/articles/nature08460#Sec3 (search "Single Sample" GSEA).

    :param pd.DataFrame data: pandas.DataFrame with the quantified features (i.e. subject x proteins)
    :param str annotation: pandas.DataFrame with the annotation to be used in the enrichment
        (i.e. CKG pathway annotation file)
    :param list[str] set_index: column/s to be used as index. Enrichment will be calculated
        for these values (i.e ["subject"] will return subjects x pathways matrix of
        enrichment scores)
    :param str annotation_col: name of the column containing annotation terms.
    :param str identifier_col: name of column containing dependent variables identifiers.
    :param str out_dir: directory path where results will be stored
        (default None, tmp folder is used)
    :param int min_size: minimum number of features (i.e. proteins) in enriched terms
        (i.e. pathways)
    :param int max_size: maximum number of features (i.e. proteins) in enriched terms
        (i.e. pathways)
    :param bool scale: whether or not to scale the data
    :param int permutations: number of permutations used in the ssgsea analysis
    :return: pandas.DataFrame containing unnormalized enrichment scores (`ES`) for each sample,
        and normalized enrichment scores (`NES`) with the enriched `Term` and sample `Name`.
    :rtype: pandas.DataFrame

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
    df = data.copy()
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Comine columns to create a unique name for each set (?)
    name = []
    if set_index is None:
        index = data.index.to_frame()
        set_index = index.columns.tolist()
    else:
        index = data[set_index]
        df = df.drop(set_index, axis=1)

    for _, row in index.iterrows():
        name.append(
            "_".join(row[set_index].tolist())
        )  # this assumes strings as identifiers

    df["Name"] = name
    index.index = name
    df = df.set_index("Name").transpose()

    if not annotation_col in annotation:
        raise ValueError(
            f"Missing Annotation Column: {annotation_col} as specified by `annotation_col`"
        )

    if not identifier_col in annotation:
        raise ValueError(
            f"Missing Identifier Column: {identifier_col} as specified by `identifier_col`"
        )

    grouped_annotations = (
        annotation.groupby(annotation_col)[identifier_col].apply(list).reset_index()
    )
    fid = uuid.uuid4()
    file_path = os.path.join(outdir, str(fid) + ".gmt")
    with open(file_path, "w", encoding="utf8") as out:
        for _, row in grouped_annotations.iterrows():
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
    result = pd.DataFrame(enrichment.res2d).set_index("Name")
    # potentially return wide format in separate format
    # result = {"es": enrichment_es, "nes": enrichment_nes}
    return result
