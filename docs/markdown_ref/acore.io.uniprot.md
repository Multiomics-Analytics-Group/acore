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

### filter_annotations(annotations: DataFrame, keywords: [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)], case_sensitive: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → DataFrame

Filter a long-format annotations DataFrame to rows whose `annotation`
column contains **any** of the given keywords.

* **Parameters:**
  * **annotations** (*pd.DataFrame*) – Long-format DataFrame with at least an `annotation` column (as
    returned by `query_uniprot()`).
  * **keywords** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *[*[*str*](https://docs.python.org/3/library/stdtypes.html#str) *]*) – Keywords to search for.  A row is kept when at least one keyword
    appears in the annotation string.  Empty strings are ignored.
  * **case_sensitive** ([*bool*](https://docs.python.org/3/library/functions.html#bool) *,* *optional*) – Whether the search is case-sensitive.  Default `False`.
* **Returns:**
  Filtered copy of *annotations* where at least one keyword matched.
  Returns an empty DataFrame with the same columns when no keywords
  are provided or none match.
* **Return type:**
  pd.DataFrame

### Examples

```pycon
>>> import pandas as pd
>>> df = pd.DataFrame({
...     "identifier": ["P1", "P1", "P2"],
...     "source": ["go_p", "go_p", "go_p"],
...     "annotation": [
...         "cell apoptosis [GO:0042981]",
...         "cell division",
...         "mitochondria",
...     ],
... })
>>> result = filter_annotations(df, keywords=["apoptosis", "mitochondr"])
>>> list(result["annotation"])
['cell apoptosis [GO:0042981]', 'mitochondria']
>>> filter_annotations(df, keywords=[]).shape[0]
0
```

## Submodules

## acore.io.uniprot.filter module

Filtering annotations by keywords.

from acore.io.uniprot import fetch_annotations, process_annotations, filter_annotations
fields = “accession,go_p,go_c,go_f”
uniprot_ids = [“P12345”, “Q67890”, “P05067”, “A1B2C3”]
df = fetch_annotations(uniprot_ids, fields=fields)
df = process_annotations(df, fields=fields)
# first keyword does not exist
filtered = filter_annotations(df, keywords=[“apoptosis”, “mitochondr”, “synapse”])

### filter_annotations(annotations: DataFrame, keywords: [list](https://docs.python.org/3/library/stdtypes.html#list)[[str](https://docs.python.org/3/library/stdtypes.html#str)], case_sensitive: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → DataFrame

Filter a long-format annotations DataFrame to rows whose `annotation`
column contains **any** of the given keywords.

* **Parameters:**
  * **annotations** (*pd.DataFrame*) – Long-format DataFrame with at least an `annotation` column (as
    returned by `query_uniprot()`).
  * **keywords** ([*list*](https://docs.python.org/3/library/stdtypes.html#list) *[*[*str*](https://docs.python.org/3/library/stdtypes.html#str) *]*) – Keywords to search for.  A row is kept when at least one keyword
    appears in the annotation string.  Empty strings are ignored.
  * **case_sensitive** ([*bool*](https://docs.python.org/3/library/functions.html#bool) *,* *optional*) – Whether the search is case-sensitive.  Default `False`.
* **Returns:**
  Filtered copy of *annotations* where at least one keyword matched.
  Returns an empty DataFrame with the same columns when no keywords
  are provided or none match.
* **Return type:**
  pd.DataFrame

### Examples

```pycon
>>> import pandas as pd
>>> df = pd.DataFrame({
...     "identifier": ["P1", "P1", "P2"],
...     "source": ["go_p", "go_p", "go_p"],
...     "annotation": [
...         "cell apoptosis [GO:0042981]",
...         "cell division",
...         "mitochondria",
...     ],
... })
>>> result = filter_annotations(df, keywords=["apoptosis", "mitochondr"])
>>> list(result["annotation"])
['cell apoptosis [GO:0042981]', 'mitochondria']
>>> filter_annotations(df, keywords=[]).shape[0]
0
```

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
