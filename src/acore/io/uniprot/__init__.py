"""Uniprot API user functions for fetching annotations for UniProt IDs and providing
the results as a pandas.DataFrame."""

import pandas as pd

from .uniprot import (
    check_id_mapping_results_ready,
    get_id_mapping_results_link,
    get_id_mapping_results_search,
    submit_id_mapping,
)


# function for outside usage
def fetch_annotations(
    ids: pd.Index | list,
    fields: str = "accession,go_p,go_c,go_f",
) -> pd.DataFrame:
    """Fetch annotations for UniProt IDs. Combines several calls to the API of UniProt's
    knowledgebase (KB).

    Parameters
    ----------
    ids : pd.Index | list
        Iterable of UniProt IDs. Fetches annotations as speecified by the specified fields.
    fields : str, optional
        Fields to fetch, by default "accession,go_p,go_c. See for availble fields:
        https://www.uniprot.org/help/return_fields

    Returns
    -------
    pd.DataFrame
        DataFrame with annotations of the UniProt IDs.
    """
    job_id = submit_id_mapping(from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=ids)
    # tsv used here to return a DataFrame. Maybe other formats are availale at some points
    _format = "tsv"
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        # add fields to the link to get more information
        # From and Entry (accession) are the same for UniProt IDs.
        results = get_id_mapping_results_search(
            link + f"?fields={fields}&format={_format}"
        )
    header = results.pop(0).split("\t")
    results = [line.split("\t") for line in results]
    df = pd.DataFrame(results, columns=header)
    return df


def process_annotations(annotations: pd.DataFrame, fields: str) -> pd.DataFrame:
    """Process annotations fetched from UniProt API.

    Parameters
    ----------
    annotations : pd.DataFrame
        DataFrame with annotations fetched from UniProt API.
    fields : str
        Fields that were fetched from the API. Comma-separated string. Fields
        needs to match number of columns in annotations.

    Returns
    -------
    pd.DataFrame
        Processed DataFrame with annotations in long-format.
    """
    d_fields_to_col = {
        k: v for k, v in zip(fields.split(","), annotations.columns[1:], strict=True)
    }

    # expand go terms
    to_expand = list()
    for field in d_fields_to_col:
        if "go_" in field:
            col = d_fields_to_col[field]
            annotations[col] = annotations[col].str.split(";")
            to_expand.append(col)
    for col in to_expand:
        # this is a bit wastefull. Processing to stack format should be done here.
        annotations = annotations.explode(col, ignore_index=True)
    # process other than go term columns
    annotations = (
        annotations.set_index("From")
        .rename_axis("identifier")
        # .drop("Entry", axis=1)
        .rename_axis("source", axis=1)
        .stack()
        .to_frame("annotation")
        .reset_index()
        .drop_duplicates(ignore_index=True)
        .replace("", pd.NA)
        .dropna()
    )
    return annotations
