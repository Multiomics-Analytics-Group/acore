"""investigate_smiles.py
This script takes a SMILES as input and extracts properties from
PubChem.

Inputs:
    - compound SMILES as string

Outputs:
    - properties: compound properties
    - data: all data from pubchem entry

Prints:
    - properties

Functions:
    - query_pubchem(smiles)
    - print_properties(properties)

"""

# %%
import time
import requests

BASE_PUBCHEM_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def inchikey_to_cid(inchikey: str) -> int | None:
    """Look up PubChem CID from an InChIKey."""
    url = f"{BASE_PUBCHEM_API}/compound/inchikey/{inchikey}/cids/JSON"
    r = requests.get(url)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    data = r.json()
    cids = data.get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None


def cid_to_kegg(cid: int) -> list[str]:
    """
    Get KEGG compound IDs for a given PubChem CID.
    Uses the KEGG REST conv endpoint, which maps pubchem SIDs → KEGG IDs.
    Then cross-references via PubChem's xrefs to get the SIDs first.
    """
    # Step 1: Get SIDs linked to this CID (KEGG deposits as SIDs)
    url = f"{BASE_PUBCHEM_API}/compound/cid/{cid}/sids/JSON"
    r = requests.get(url)
    if r.status_code == 404:
        return []
    r.raise_for_status()
    sids = (
        r.json().get("InformationList", {}).get("Information", [{}])[0].get("SID", [])
    )
    print(f"  Found SIDs for CID {cid}: {sids}")
    kegg_ids = []
    for sid in sids:

        kegg_url = f"https://rest.kegg.jp/conv/compound/pubchem:{sid}"
        kr = requests.get(kegg_url)
        if kr.status_code == 200 and kr.text.strip():
            for line in kr.text.strip().splitlines():
                parts = line.split("\t")
                if len(parts) == 2:
                    kegg_ids.append(parts[1])  # e.g. "cpd:C00031"
                if len(kegg_ids) > 0:
                    return kegg_ids
        time.sleep(0.2)  # be polite to KEGG

    return list(set(kegg_ids))


if __name__ == "__main__":
    # %%
    inchikey = "WQZGKKKJIJ
