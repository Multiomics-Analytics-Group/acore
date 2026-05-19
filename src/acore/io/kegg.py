# %%
# region: Imports and constants
import re
from collections.abc import Iterable
from typing import Iterable
from urllib import error, request

import pandas as pd
import requests

KEGG_API_BASE_URL = "https://rest.kegg.jp"
_KO_TERM_PATTERN = re.compile(r"^(?:ko:)?(K\d{5})$", re.IGNORECASE)
MAX_KEGG_BATCH_SIZE = 10


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
    for i in range(0, len(gene_ids), MAX_KEGG_BATCH_SIZE):
        batch = "+".join(gene_ids[i : i + MAX_KEGG_BATCH_SIZE])
        r = requests.get(f"https://rest.kegg.jp/link/{target_db}/{batch}", timeout=30)
        r.raise_for_status()
        results.append(r.text)
    return "".join(results)


def fetch_kegg_ko_descriptions(
    ko_terms: Iterable[str],
    timeout: float = 30.0,
) -> pd.DataFrame:
    """Fetch common descriptions for KEGG KO terms.

    Parameters
    ----------
    ko_terms : Iterable[str]
        KEGG KO identifiers such as ``K03007`` or ``ko:K03007``.
    timeout : float, optional
        Timeout in seconds for each KEGG API request.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns ``ko_term``, ``symbol`` and
        ``common_description``.

    Notes
    -----
    The KEGG API accepts up to 10 entry identifiers per request. This helper
    batches larger inputs automatically.
    """
    normalized_terms = _unique_in_order(_normalize_ko_term(term) for term in ko_terms)
    if not normalized_terms:
        return pd.DataFrame(columns=["ko_term", "symbol", "common_description"]).astype(
            str
        )

    rows: list[dict[str, str]] = []
    for batch in _batched(normalized_terms, MAX_KEGG_BATCH_SIZE):
        response_text = _fetch_kegg_batch(batch, timeout=timeout)
        rows.extend(_parse_kegg_ko_entries(response_text))

    return pd.DataFrame(rows, columns=["ko_term", "symbol", "common_description"])


def _fetch_kegg_batch(ko_terms: list[str], timeout: float) -> str:
    batch_arg = "+".join(ko_terms)
    url = f"{KEGG_API_BASE_URL}/get/{batch_arg}"
    try:
        with request.urlopen(url, timeout=timeout) as response:
            return response.read().decode("utf-8")
    except error.URLError as exc:
        raise RuntimeError(f"Failed to fetch KEGG KO data from {url}") from exc


def _parse_kegg_ko_entries(response_text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for entry_text in response_text.split("///"):
        stripped_entry = entry_text.strip()
        if not stripped_entry:
            continue
        fields = _parse_kegg_flat_file_entry(stripped_entry)
        entry_id = fields.get("ENTRY", "").split()[0]
        if not entry_id:
            continue
        rows.append(
            {
                "ko_term": f"ko:{entry_id}",
                "symbol": fields.get("SYMBOL", ""),
                "common_description": fields.get("NAME", ""),
            }
        )
    return rows


def _parse_kegg_flat_file_entry(entry_text: str) -> dict[str, str]:
    parsed_fields: dict[str, list[str]] = {}
    current_field = ""

    for line in entry_text.splitlines():
        field_name = line[:12].strip()
        field_value = line[12:].strip()
        if field_name:
            current_field = field_name
            parsed_fields.setdefault(current_field, []).append(field_value)
            continue
        if current_field:
            parsed_fields[current_field].append(field_value)

    return {
        field_name: " ".join(values).strip()
        for field_name, values in parsed_fields.items()
    }


def _normalize_ko_term(ko_term: str) -> str:
    stripped_term = ko_term.strip()
    match = _KO_TERM_PATTERN.fullmatch(stripped_term)
    if not match:
        raise ValueError(
            f"Invalid KEGG KO term {ko_term!r}. Expected values like 'K03007' or "
            "'ko:K03007'."
        )
    return f"ko:{match.group(1).upper()}"


def _unique_in_order(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(values))


def _batched(values: list[str], batch_size: int) -> Iterable[list[str]]:
    for index in range(0, len(values), batch_size):
        yield values[index : index + batch_size]


# endregion

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
    # https://rest.kegg.jp/list/path:hsa00010+path:hsa00500
    # curl -fsSL https://rest.kegg.jp/link/pathway/yli:2912002
    # curl -fsSL https://rest.kegg.jp/get/path:yli00592
    # %%
    df = fetch_kegg_ko_descriptions(["ko:K03007", "K02143", "ko:K00844"])
    df
# %%
