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
