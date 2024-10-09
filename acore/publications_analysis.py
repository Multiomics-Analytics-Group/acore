from urllib import error

import pandas as pd
from Bio import Entrez, Medline

# TODO: This should probably be changed to the email of the person installing ckg?
Entrez.email = "kg@dtu.dk"


def getMedlineAbstracts(idList):
    fields = {
        "TI": "title",
        "AU": "authors",
        "JT": "journal",
        "DP": "date",
        "MH": "keywords",
        "AB": "abstract",
        "PMID": "PMID",
    }
    pubmedUrl = "https://www.ncbi.nlm.nih.gov/pubmed/"
    abstracts = pd.DataFrame()
    try:
        handle = Entrez.efetch(
            db="pubmed", id=idList, rettype="medline", retmode="json"
        )
        records = Medline.parse(handle)
        results = []
        for record in records:
            aux = {}
            for field in fields:
                if field in record:
                    aux[fields[field]] = record[field]

            if "PMID" in aux:
                aux["url"] = pubmedUrl + aux["PMID"]
            else:
                aux["url"] = ""
            results.append(aux)

        abstracts = pd.DataFrame.from_dict(results)
    except error.URLError as e:
        print("URLError: Request to Bio.Entrez failed. Error: {}".format(e))
    except error.HTTPError as e:
        print("HTTPError: Request to Bio.Entrez failed. Error: {}".format(e))
    except Exception as e:
        print("Request to Bio.Entrez failed. Error: {}".format(e))

    return abstracts


def get_publications_abstracts(
    data,
    publication_col="publication",
    join_by=["publication", "Proteins", "Diseases"],
    index="PMID",
):
    """
    Accesses NCBI PubMed over the WWW and retrieves the abstracts corresponding
    to a list of one or more PubMed IDs.

    :param data: pandas dataframe of diseases and publications linked to a list of
                proteins (columns: 'Diseases', 'Proteins', 'linkout' and 'publication').
    :param str publication_col: column label containing PubMed ids.
    :param list join_by: column labels to be kept from the input dataframe.
    :param str index: column label containing PubMed ids from the NCBI retrieved data.
    :return: Pandas dataframe with publication information and columns 'PMID', 'abstract',
     'authors', 'date', 'journal', 'keywords', 'title', 'url', 'Proteins' and 'Diseases'.

    Example::

        result = get_publications_abstracts(data,
                    publication_col='publication',
                    join_by=['publication','Proteins','Diseases'],
                    index='PMID')
    """
    abstracts = pd.DataFrame()
    if not data.empty:
        abstracts = getMedlineAbstracts(
            list(data.reset_index()[publication_col].unique())
        )
        if not abstracts.empty:
            abstracts = abstracts.set_index(index)
            abstracts = abstracts.join(
                data.reset_index()[join_by].set_index(publication_col)
            )
            abstracts.index.name = index
            abstracts = abstracts.reset_index()
    return abstracts
