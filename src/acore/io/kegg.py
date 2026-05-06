from typing import Iterable

import requests


def link_kegg_batch(target_db: str, gene_ids: Iterable[str]) -> str:
    """Fetch from KEGG in batches informations.

    Docs: https://www.kegg.jp/kegg/rest/keggapi.html

    Parameters
    ----------
    target_db : str
        Target endpoint of KEGG API, e.g. "ko" or "pathway".
    gene_ids : list of str
        List of KEGG gene IDs to query.

    Returns
    -------
    _type_
        _description_
    """

    results = []
    # maximum of 10 gene IDs per request, as per KEGG API documentation for GET requests
    # ! To Do: see if this is true also for link requests.
    for i in range(0, len(gene_ids), 10):
        batch = "+".join(gene_ids[i : i + 10])
        r = requests.get(f"https://rest.kegg.jp/link/{target_db}/{batch}", timeout=30)
        r.raise_for_status()
        results.append(r.text)
    return "".join(results)


if __name__ == "__main__":
    # Example KEGG gene IDs for testing (fetched from UniProt)
    kegg_genes = [
        "yli:7009397",
        "yli:7009411",
        "yli:2910294",
        "yli:2912002",
        "hsa:3099",
    ]
    results = link_kegg_batch("ko", kegg_genes)
    print(results)
    pathways = link_kegg_batch("pathway", kegg_genes)
    print(pathways)
