# acore.io.uniprot package

Uniprot API user functions for fetching annotations for UniProt IDs and providing
the results as a pandas.DataFrame.

### fetch_annotations(ids: Index | [list](https://docs.python.org/3/library/stdtypes.html#list), fields: [str](https://docs.python.org/3/library/stdtypes.html#str) = 'accession,go_p,go_c,go_f') → DataFrame

Fetch annotations for UniProt IDs. Combines several calls to the API of UniProt’s
knowledgebase (KB).

* **Parameters:**
  * **ids** (*pd.Index* *|* [*list*](https://docs.python.org/3/library/stdtypes.html#list)) – Iterable of UniProt IDs. Fetches annotations as speecified by the specified fields.
  * **fields** ([*str*](https://docs.python.org/3/library/stdtypes.html#str) *,* *optional*) – Fields to fetch, by default “accession,go_p,go_c. See for availble fields:
    [https://www.uniprot.org/help/return_fields](https://www.uniprot.org/help/return_fields)
* **Returns:**
  DataFrame with annotations of the UniProt IDs.
* **Return type:**
  pd.DataFrame

### process_annotations(annotations: DataFrame, fields: [str](https://docs.python.org/3/library/stdtypes.html#str)) → DataFrame

Process annotations fetched from UniProt API.

* **Parameters:**
  * **annotations** (*pd.DataFrame*) – DataFrame with annotations fetched from UniProt API.
  * **fields** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Fields that were fetched from the API. Comma-separated string. Fields
    needs to match number of columns in annotations.
* **Returns:**
  Processed DataFrame with annotations in long-format.
* **Return type:**
  pd.DataFrame

## Submodules

## acore.io.uniprot.uniprot module

Uniprot ID mapping using Python.

Source: [https://www.uniprot.org/help/id_mapping](https://www.uniprot.org/help/id_mapping)

### check_response(response)

### submit_id_mapping(from_db, to_db, ids)

### get_next_link(headers)

### check_id_mapping_results_ready(job_id)

### get_batch(batch_response, file_format, compressed)

### combine_batches(all_results, batch_results, file_format)

### get_id_mapping_results_link(job_id)

### decode_results(response, file_format, compressed)

### get_xml_namespace(element)

### merge_xml_results(xml_results)

### print_progress_batches(batch_index, size, total)

### get_id_mapping_results_search(url)

### get_id_mapping_results_stream(url)
