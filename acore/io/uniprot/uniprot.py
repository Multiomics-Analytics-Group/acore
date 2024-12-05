"""Uniprot ID mapping using Python.

Source: https://www.uniprot.org/help/id_mapping
"""

import json
import re
import time
import zlib
from urllib.parse import parse_qs, urlencode, urlparse
from xml.etree import ElementTree

import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry

POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"


retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] in ("NEW", "RUNNING"):
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")


def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results


def get_id_mapping_results_stream(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/results/stream/")
    request = session.get(url)
    check_response(request)
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)




if __name__ == "__main__":
    # id mapping is used to create a link to a query (you can see the json in the browser)
    # UniProtKB is the knowleadgebase integrating all kind of other databases
    import pandas as pd

    job_id = submit_id_mapping(
        from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=["P05067", "P12345"]
    )

    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        # add fields to the link to get more information
        # From and Entry (accession) are the same.
        results = get_id_mapping_results_search(
            link + "?fields=accession,go_p,go_c,go_f&format=tsv"
        )
        # see the available fields you can add
        # https://www.uniprot.org/help/return_fields
        # and the available formats: json and tsv (for most endpoints)
        # https://www.uniprot.org/help/api_queries#tips
    print(results)
    header = results.pop(0).split("\t")
    results = [line.split("\t") for line in results]
    df = pd.DataFrame(results, columns=header)
    # result from:
    # df.to_dict(orient="records")
    records = [
        {
            "From": "P05067",
            "Gene Ontology (biological process)": (
                "adult locomotory behavior [GO:0008344]; amyloid fibril formation"
                " [GO:1990000]; astrocyte activation [GO:0048143]; astrocyte activation"
                " involved in immune response [GO:0002265]; axo-dendritic transport"
                " [GO:0008088]; axon midline choice point recognition [GO:0016199];"
                " axonogenesis [GO:0007409]; cell adhesion [GO:0007155]; cellular"
                " response to amyloid-beta [GO:1904646]; central nervous system"
                " development [GO:0007417]; cholesterol metabolic process [GO:0008203];"
                " cognition [GO:0050890]; collateral sprouting in absence of injury"
                " [GO:0048669]; cytosolic mRNA polyadenylation [GO:0180011]; dendrite"
                " development [GO:0016358]; endocytosis [GO:0006897]; extracellular"
                " matrix organization [GO:0030198]; forebrain development [GO:0030900];"
                " G2/M transition of mitotic cell cycle [GO:0000086]; intracellular"
                " copper ion homeostasis [GO:0006878]; ionotropic glutamate receptor"
                " signaling pathway [GO:0035235]; learning [GO:0007612]; learning or"
                " memory [GO:0007611]; locomotory behavior [GO:0007626]; mating"
                " behavior [GO:0007617]; microglia development [GO:0014005]; microglial"
                " cell activation [GO:0001774]; modulation of excitatory postsynaptic"
                " potential [GO:0098815]; negative regulation of cell population"
                " proliferation [GO:0008285]; negative regulation of gene expression"
                " [GO:0010629]; negative regulation of long-term synaptic potentiation"
                " [GO:1900272]; negative regulation of neuron differentiation"
                " [GO:0045665]; neuromuscular process controlling balance [GO:0050885];"
                " neuron apoptotic process [GO:0051402]; neuron cellular homeostasis"
                " [GO:0070050]; neuron projection development [GO:0031175]; neuron"
                " projection maintenance [GO:1990535]; neuron remodeling [GO:0016322];"
                " NMDA selective glutamate receptor signaling pathway [GO:0098989];"
                " Notch signaling pathway [GO:0007219]; positive regulation of amyloid"
                " fibril formation [GO:1905908]; positive regulation of"
                " calcium-mediated signaling [GO:0050850]; positive regulation of"
                " chemokine production [GO:0032722]; positive regulation of ERK1 and"
                " ERK2 cascade [GO:0070374]; positive regulation of G2/M transition of"
                " mitotic cell cycle [GO:0010971]; positive regulation of gene"
                " expression [GO:0010628]; positive regulation of glycolytic process"
                " [GO:0045821]; positive regulation of inflammatory response"
                " [GO:0050729]; positive regulation of interleukin-1 beta production"
                " [GO:0032731]; positive regulation of interleukin-6 production"
                " [GO:0032755]; positive regulation of JNK cascade [GO:0046330];"
                " positive regulation of long-term synaptic potentiation [GO:1900273];"
                " positive regulation of mitotic cell cycle [GO:0045931]; positive"
                " regulation of non-canonical NF-kappaB signal transduction"
                " [GO:1901224]; positive regulation of peptidyl-serine phosphorylation"
                " [GO:0033138]; positive regulation of peptidyl-threonine"
                " phosphorylation [GO:0010800]; positive regulation of protein"
                " metabolic process [GO:0051247]; positive regulation of protein"
                " phosphorylation [GO:0001934]; positive regulation of T cell migration"
                " [GO:2000406]; positive regulation of transcription by RNA polymerase"
                " II [GO:0045944]; positive regulation of tumor necrosis factor"
                " production [GO:0032760]; regulation of gene expression [GO:0010468];"
                " regulation of long-term neuronal synaptic plasticity [GO:0048169];"
                " regulation of multicellular organism growth [GO:0040014]; regulation"
                " of peptidyl-tyrosine phosphorylation [GO:0050730]; regulation of"
                " presynapse assembly [GO:1905606]; regulation of spontaneous synaptic"
                " transmission [GO:0150003]; regulation of synapse structure or"
                " activity [GO:0050803]; regulation of translation [GO:0006417];"
                " regulation of Wnt signaling pathway [GO:0030111]; response to"
                " interleukin-1 [GO:0070555]; response to oxidative stress"
                " [GO:0006979]; smooth endoplasmic reticulum calcium ion homeostasis"
                " [GO:0051563]; suckling behavior [GO:0001967]; synapse organization"
                " [GO:0050808]; synaptic assembly at neuromuscular junction"
                " [GO:0051124]; visual learning [GO:0008542]"
            ),
            "Gene Ontology (cellular component)": (
                "apical part of cell [GO:0045177]; axon [GO:0030424]; cell surface"
                " [GO:0009986]; cell-cell junction [GO:0005911]; ciliary rootlet"
                " [GO:0035253]; clathrin-coated pit [GO:0005905]; COPII-coated ER to"
                " Golgi transport vesicle [GO:0030134]; cytoplasm [GO:0005737]; cytosol"
                " [GO:0005829]; dendrite [GO:0030425]; dendritic shaft [GO:0043198];"
                " dendritic spine [GO:0043197]; early endosome [GO:0005769];"
                " endoplasmic reticulum [GO:0005783]; endoplasmic reticulum lumen"
                " [GO:0005788]; endosome [GO:0005768]; endosome lumen [GO:0031904];"
                " extracellular exosome [GO:0070062]; extracellular region"
                " [GO:0005576]; extracellular space [GO:0005615]; Golgi apparatus"
                " [GO:0005794]; Golgi lumen [GO:0005796]; Golgi-associated vesicle"
                " [GO:0005798]; growth cone [GO:0030426]; membrane [GO:0016020];"
                " membrane raft [GO:0045121]; mitochondrial inner membrane"
                " [GO:0005743]; neuromuscular junction [GO:0031594]; nuclear envelope"
                " lumen [GO:0005641]; perikaryon [GO:0043204]; perinuclear region of"
                " cytoplasm [GO:0048471]; plasma membrane [GO:0005886]; platelet alpha"
                " granule lumen [GO:0031093]; presynaptic active zone [GO:0048786];"
                " receptor complex [GO:0043235]; recycling endosome [GO:0055037];"
                " smooth endoplasmic reticulum [GO:0005790]; spindle midzone"
                " [GO:0051233]; synapse [GO:0045202]; synaptic vesicle [GO:0008021];"
                " trans-Golgi network membrane [GO:0032588]"
            ),
            "Gene Ontology (molecular function)": (
                "DNA binding [GO:0003677]; enzyme binding [GO:0019899]; heparin binding"
                " [GO:0008201]; identical protein binding [GO:0042802]; protein"
                " serine/threonine kinase binding [GO:0120283]; PTB domain binding"
                " [GO:0051425]; receptor ligand activity [GO:0048018]; RNA polymerase"
                " II cis-regulatory region sequence-specific DNA binding [GO:0000978];"
                " serine-type endopeptidase inhibitor activity [GO:0004867]; signaling"
                " receptor activator activity [GO:0030546]; signaling receptor binding"
                " [GO:0005102]; transition metal ion binding [GO:0046914]"
            ),
        },
        {
            "From": "P12345",
            "Gene Ontology (biological process)": (
                "2-oxoglutarate metabolic process [GO:0006103]; aspartate catabolic"
                " process [GO:0006533]; aspartate metabolic process [GO:0006531];"
                " glutamate metabolic process [GO:0006536]; lipid transport"
                " [GO:0006869]; protein folding [GO:0006457]"
            ),
            "Gene Ontology (cellular component)": (
                "mitochondrial matrix [GO:0005759]; mitochondrion [GO:0005739]; plasma"
                " membrane [GO:0005886]"
            ),
            "Gene Ontology (molecular function)": (
                "kynurenine-oxoglutarate transaminase activity [GO:0016212];"
                " L-aspartate:2-oxoglutarate aminotransferase activity [GO:0004069];"
                " pyridoxal phosphate binding [GO:0030170]"
            ),
        },
    ]
