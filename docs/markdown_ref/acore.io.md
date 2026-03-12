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

### unrar(filepath, to)

Decompress RAR file
:param str filepath: path to rar file
:param str to: where to extract all files

### download_file(url: [str](https://docs.python.org/3/library/stdtypes.html#str), local_filename: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [None](https://docs.python.org/3/library/constants.html#None)

Download a file from the internet.

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

## acore.io.ftp module

### download_from_ftp(ftp_url: [str](https://docs.python.org/3/library/stdtypes.html#str), user: [str](https://docs.python.org/3/library/stdtypes.html#str), password: [str](https://docs.python.org/3/library/stdtypes.html#str), to: [str](https://docs.python.org/3/library/stdtypes.html#str), file_name) → [str](https://docs.python.org/3/library/stdtypes.html#str)

Download a file from an FTP server.

## acore.io.http module

Download files from the internet.

### download_file(url: [str](https://docs.python.org/3/library/stdtypes.html#str), local_filename: [str](https://docs.python.org/3/library/stdtypes.html#str)) → [None](https://docs.python.org/3/library/constants.html#None)

Download a file from the internet.

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
