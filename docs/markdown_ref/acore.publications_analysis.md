# acore.publications_analysis package

### getMedlineAbstracts(idList)

### get_publications_abstracts(data, publication_col='publication', join_by=['publication', 'Proteins', 'Diseases'], index='PMID')

Accesses NCBI PubMed over the WWW and retrieves the abstracts corresponding
to a list of one or more PubMed IDs.

* **Parameters:**
  * **data** – pandas dataframe of diseases and publications linked to a list of
    proteins (columns: ‘Diseases’, ‘Proteins’, ‘linkout’ and ‘publication’).
  * **publication_col** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing PubMed ids.
  * **join_by** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)) – column labels to be kept from the input dataframe.
  * **index** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – column label containing PubMed ids from the NCBI retrieved data.
* **Returns:**
  Pandas dataframe with publication information and columns ‘PMID’, ‘abstract’,
  ‘authors’, ‘date’, ‘journal’, ‘keywords’, ‘title’, ‘url’, ‘Proteins’ and ‘Diseases’.

Example:

```default
result = get_publications_abstracts(data,
            publication_col='publication',
            join_by=['publication','Proteins','Diseases'],
            index='PMID')
```
