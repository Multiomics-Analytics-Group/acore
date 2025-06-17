"""Enrichment Analysis Module. Contains different functions to perform enrichment
analysis.

Most things in this module are covered in https://www.youtube.com/watch?v=2NC1QOXmc5o
by Lars Juhl Jensen.
"""

from __future__ import annotations

import logging
import os
import re
import uuid

import gseapy as gp
import pandas as pd

from acore.enrichment_analysis.annotate import annotate_features
from acore.enrichment_analysis.statistical_tests.fisher import run_fisher
from acore.enrichment_analysis.statistical_tests.kolmogorov_smirnov import (
    run_kolmogorov_smirnov,
)
from acore.multiple_testing import apply_pvalue_correction

logger = logging.getLogger(__name__)

TYPE_COLS_MSG = """
columns: 'terms', 'identifiers', 'foreground',
    'background', foreground_pop, background_pop, 'pvalue', 'padj' and 'rejected'.
"""

__all__ = [
    "run_site_regulation_enrichment",
    "run_up_down_regulation_enrichment",
    "run_fisher",
    "run_kolmogorov_smirnov",
]


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
        regulation_data=regulation_data,
        annotation=annotation,
        identifier=identifier,
        group_col=groups,
        annotation_col=annotation_col,
        rejected_col=rejected_col,
        method=method,
        correction=correction,
    )

    return result


def run_up_down_regulation_enrichment(
    regulation_data: pd.DataFrame,
    annotation: pd.DataFrame,
    identifier: str = "identifier",
    groups: list[str] = ("group1", "group2"),
    annotation_col: str = "annotation",
    # rejected_col: str = "rejected", # could be passed
    pval_col: str = "pval",
    group_col: str = "group",
    log2fc_col: str = "log2FC",
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

        df["up_pairwise_regulation"] = (df[pval_col] <= correction_alpha) & (
            df[log2fc_col] >= lfc_cutoff
        )
        df["down_pairwise_regulation"] = (df[pval_col] <= correction_alpha) & (
            df[log2fc_col] <= -lfc_cutoff
        )
        comparison_tag = str(g1) + "~" + str(g2)

        if not regulation_data[identifier].is_unique:
            logger.warning(
                "Column '%s' in regulation_data contains duplicated values for comparison %s.",
                identifier,
                comparison_tag,
            )

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
        It should be present both in `regulation_data` and `annotation`. In `regulation_data`
        it should be unique, while in `annotation` it can contain duplicates as one
        identifier can be part of multiple pathways.
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
    if regulation_data[rejected_col].isna().any():
        raise ValueError(f"Rejected column '{rejected_col}' contains missing values.")
    mask_rejected = regulation_data[rejected_col].astype(bool)
    if not regulation_data[identifier].is_unique:
        raise ValueError(f"Column '{identifier}' in regulation_data has to be unique.")
    foreground_list = regulation_data.loc[mask_rejected, identifier]
    background_list = regulation_data.loc[~mask_rejected, identifier]
    foreground_pop = len(foreground_list)
    background_pop = len(regulation_data[identifier])
    # needs to allow for missing annotations
    # ! this step needs unique identifiers in the regulation_data
    # group_col contains either 'foreground', 'background' or NA
    annotation[group_col] = annotate_features(
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
    :param int min_detected_in_set: minimum number of features in the foreground
    :param str annotation_col: name of the column containing annotation terms.
    :param str group_col: name of column containing the group identifiers,
                          e.g. specifying belonging to 'foreground' or 'background'.
    :param str identifier_col: name of column containing dependent variables identifiers.
    :param str method: method used to compute enrichment
                       (only 'fisher' is supported currently).
    :param str correction: method to be used for multiple-testing correction.
    :param float correction_alpha: adjusted p-value cutoff to define significance.
    :return: pandas.DataFrame with columns: annotation terms, features,
        number of foreground/background features in each term,
        p-values and corrected p-values
        Columns are: 'terms', 'identifiers',
        'foreground', 'background', 'foreground_pop', 'background_pop',
        'pvalue', 'padj' and 'rejected'.

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
    mask_in_foreground = countsdf[group_col] == foreground_id
    terms_in_foreground = countsdf.loc[mask_in_foreground, annotation_col].unique()
    for annotation in terms_in_foreground:
        counts = countsdf[countsdf[annotation_col] == annotation]
        num_foreground = int(
            counts.loc[counts[group_col] == foreground_id, "count"].squeeze()
        )
        num_background = 0  # initialize to 0 in case all features are foreground
        num_background = counts.loc[
            counts[group_col] == background_id, "count"
        ].squeeze()
        if isinstance(num_background, pd.Series) and num_background.empty:
            # if no value is found in counts, an empty series is returned
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
                "foreground_pop": foreground_pop,  # no. of foreground features, constant
                "background_pop": background_pop,  # no. of included features, constant
                "pvalue": pvalues,
                "padj": padj,
                "rejected": rejected.astype(bool),
            }
        )
        result = result.sort_values(by="padj", ascending=True)
    else:
        logger.warning(
            "No significant enrichment found with the given parameters. "
            "Returning an empty DataFrame."
        )
        # ToDo: Should we return an empty DataFrame with the expected columns?
        result = pd.DataFrame()

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

    :param pd.DataFrame data: pandas.DataFrame with the quantified features
                              (i.e. subject x proteins)
    :param str annotation: pandas.DataFrame with the annotation to be used in the
                            enrichment (i.e. CKG pathway annotation file)
    :param list[str] set_index: column/s to be used as index. Enrichment will be
        calculated for these values (i.e ["subject"] will return subjects x pathways
        matrix of enrichment scores)
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
    :return: pandas.DataFrame containing unnormalized enrichment scores (`ES`)
             for each sample, and normalized enrichment scores (`NES`)
             with the enriched `Term` and sample `Name`.
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
            f"Missing Annotation Column: {annotation_col}"
            " as specified by `annotation_col`"
        )

    if not identifier_col in annotation:
        raise ValueError(
            f"Missing Identifier Column: {identifier_col}"
            " as specified by `identifier_col`"
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
