# acore.io package

### download_PRIDE_data(pxd_id, file_name, to='.', user='', password='', date_field='publicationDate') → [dict](https://docs.python.org/3/library/stdtypes.html#dict)

This function downloads a project file from the PRIDE repository. To see more of the
pride API, have a look at
[https://www.ebi.ac.uk/pride/ws/archive/v3/webjars/swagger-ui/index.html](https://www.ebi.ac.uk/pride/ws/archive/v3/webjars/swagger-ui/index.html)
or EBI’s commandline tool pridepy
[https://github.com/PRIDE-Archive/pridepy](https://github.com/PRIDE-Archive/pridepy)

* **Parameters:**
  * **pxd_id** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – PRIDE project identifier (id. PXD013599).
  * **file_name** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the file to dowload
  * **to** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – local directory where the file should be downloaded
  * **user** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – username to access biomedical database server if required.
  * **password** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – password to access biomedical database server if required.
  * **date_field** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – projects deposited in PRIDE are search based on date, either
    submissionData or publicationDate (default)

### download_file(url: [str](https://docs.python.org/3/library/stdtypes.html#str), local_filename: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [None](https://docs.python.org/3/library/constants.html#None)

Download a file from the internet.

### unrar(filepath, to)

Decompress RAR file
:param str filepath: path to rar file
:param str to: where to extract all files

## Subpackages

* [acore.io.uniprot package](acore.io.uniprot.md)
  * [`fetch_annotations()`](acore.io.uniprot.md#acore.io.uniprot.fetch_annotations)
  * [`process_annotations()`](acore.io.uniprot.md#acore.io.uniprot.process_annotations)
  * [Submodules](acore.io.uniprot.md#submodules)
  * [acore.io.uniprot.uniprot module](acore.io.uniprot.md#module-acore.io.uniprot.uniprot)
    * [`check_response()`](acore.io.uniprot.md#acore.io.uniprot.uniprot.check_response)
    * [`submit_id_mapping()`](acore.io.uniprot.md#acore.io.uniprot.uniprot.submit_id_mapping)
    * [`get_next_link()`](acore.io.uniprot.md#acore.io.uniprot.uniprot.get_next_link)
    * [`check_id_mapping_results_ready()`](acore.io.uniprot.md#acore.io.uniprot.uniprot.check_id_mapping_results_ready)
    * [`get_batch()`](acore.io.uniprot.md#acore.io.uniprot.uniprot.get_batch)
    * [`combine_batches()`](acore.io.uniprot.md#acore.io.uniprot.uniprot.combine_batches)
    * [`get_id_mapping_results_link()`](acore.io.uniprot.md#acore.io.uniprot.uniprot.get_id_mapping_results_link)
    * [`decode_results()`](acore.io.uniprot.md#acore.io.uniprot.uniprot.decode_results)
    * [`get_xml_namespace()`](acore.io.uniprot.md#acore.io.uniprot.uniprot.get_xml_namespace)
    * [`merge_xml_results()`](acore.io.uniprot.md#acore.io.uniprot.uniprot.merge_xml_results)
    * [`print_progress_batches()`](acore.io.uniprot.md#acore.io.uniprot.uniprot.print_progress_batches)
    * [`get_id_mapping_results_search()`](acore.io.uniprot.md#acore.io.uniprot.uniprot.get_id_mapping_results_search)
    * [`get_id_mapping_results_stream()`](acore.io.uniprot.md#acore.io.uniprot.uniprot.get_id_mapping_results_stream)

## Submodules

## acore.io.download module

Download files from the internet.

### download_file(url: [str](https://docs.python.org/3/library/stdtypes.html#str), local_filename: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [None](https://docs.python.org/3/library/constants.html#None)

Download a file from the internet.

## acore.io.ftp module

### download_from_ftp(ftp_url: [str](https://docs.python.org/3/library/stdtypes.html#str), user: [str](https://docs.python.org/3/library/stdtypes.html#str), password: [str](https://docs.python.org/3/library/stdtypes.html#str), to: [str](https://docs.python.org/3/library/stdtypes.html#str), file_name) → [str](https://docs.python.org/3/library/stdtypes.html#str)

Download a file from an FTP server.

## acore.io.kegg module

### cid_to_kegg_id(pubchem_cid: [int](https://docs.python.org/3/library/functions.html#int)) → [str](https://docs.python.org/3/library/stdtypes.html#str) | [None](https://docs.python.org/3/library/constants.html#None)

Convert a single PubChem CID to a KEGG compound ID via KEGG conv API.

### fetch_kegg_ko_descriptions(ko_terms: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/library/stdtypes.html#str)], timeout: [float](https://docs.python.org/3/library/functions.html#float) = 30.0) → DataFrame

Fetch common descriptions for KEGG KO terms.

* **Parameters:**
  * **ko_terms** (*Iterable* *[*[*str*](https://docs.python.org/3/library/stdtypes.html#str) *]*) – KEGG KO identifiers such as `K03007` or `ko:K03007`.
  * **timeout** ([*float*](https://docs.python.org/3/library/functions.html#float) *,* *optional*) – Timeout in seconds for each KEGG API request.
* **Returns:**
  A DataFrame with columns `ko_term`, `symbol` and
  `common_description`.
* **Return type:**
  pd.DataFrame

### Notes

The KEGG API accepts up to 10 entry identifiers per request. This helper
batches larger inputs automatically.

### link_kegg_batch(target_db: [str](https://docs.python.org/3/library/stdtypes.html#str), gene_ids: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[str](https://docs.python.org/3/library/stdtypes.html#str)]) → [str](https://docs.python.org/3/library/stdtypes.html#str)

Fetch from KEGG in batches informations.

Docs: [https://www.kegg.jp/kegg/rest/keggapi.html](https://www.kegg.jp/kegg/rest/keggapi.html)

* **Parameters:**
  * **target_db** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Target endpoint of KEGG API, e.g. “ko” or “pathway”.
  * **gene_ids** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *of* [*str*](https://docs.python.org/3/library/stdtypes.html#str)) – List of KEGG gene IDs to query.
* **Returns:**
  A list of strings containing the fetched information.
* **Return type:**
  [str](https://docs.python.org/3/library/stdtypes.html#str)

### lookup_cid_to_kegg_id(pubchem_cid: [Iterable](https://docs.python.org/3/library/collections.abc.html#collections.abc.Iterable)[[int](https://docs.python.org/3/library/functions.html#int)]) → Series | [None](https://docs.python.org/3/library/constants.html#None)

Look up KEGG IDs for a list of PubChem CIDs using a pre-downloaded mapping file.

* **Parameters:**
  **pubchem_cid** (*Iterable* *[*[*int*](https://docs.python.org/3/library/functions.html#int) *]*) – A list of PubChem CIDs to look up.
* **Returns:**
  A Series mapping PubChem CIDs to KEGG IDs, or None if no matches are found.
* **Return type:**
  pd.Series | None

### parse_compound_pathway_mapping(raw_mapping: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)]]

Parse tab-delimited KEGG-style compound/pathway mappings into a dictionary.

### parse_kegg_name_description(raw_text: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [dict](https://docs.python.org/3/library/stdtypes.html#dict)[[str](https://docs.python.org/3/library/stdtypes.html#str), [str](https://docs.python.org/3/library/stdtypes.html#str)]]

Parse KEGG pathway entries into ENTRY -> {NAME, DESCRIPTION}.

## acore.io.pride module

Downlaod data from PRIDE database.

### download_PRIDE_data(pxd_id, file_name, to='.', user='', password='', date_field='publicationDate') → [dict](https://docs.python.org/3/library/stdtypes.html#dict)

This function downloads a project file from the PRIDE repository. To see more of the
pride API, have a look at
[https://www.ebi.ac.uk/pride/ws/archive/v3/webjars/swagger-ui/index.html](https://www.ebi.ac.uk/pride/ws/archive/v3/webjars/swagger-ui/index.html)
or EBI’s commandline tool pridepy
[https://github.com/PRIDE-Archive/pridepy](https://github.com/PRIDE-Archive/pridepy)

* **Parameters:**
  * **pxd_id** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – PRIDE project identifier (id. PXD013599).
  * **file_name** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – name of the file to dowload
  * **to** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – local directory where the file should be downloaded
  * **user** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – username to access biomedical database server if required.
  * **password** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – password to access biomedical database server if required.
  * **date_field** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – projects deposited in PRIDE are search based on date, either
    submissionData or publicationDate (default)

## acore.io.uncompress module

### unrar(filepath, to)

Decompress RAR file
:param str filepath: path to rar file
:param str to: where to extract all files
