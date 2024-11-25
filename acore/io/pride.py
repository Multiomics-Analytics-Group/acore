"""Downlaod data from PRIDE database."""

import requests

from .ftp import download_from_ftp


def download_PRIDE_data(
    pxd_id, file_name, to=".", user="", password="", date_field="publicationDate"
) -> dict:
    """
    This function downloads a project file from the PRIDE repository. To see more of the
    pride API, have a look at
    https://www.ebi.ac.uk/pride/ws/archive/v3/webjars/swagger-ui/index.html
    or EBI's commandline tool pridepy
    https://github.com/PRIDE-Archive/pridepy

    :param str pxd_id: PRIDE project identifier (id. PXD013599).
    :param str file_name: name of the file to dowload
    :param str to: local directory where the file should be downloaded
    :param str user: username to access biomedical database server if required.
    :param str password: password to access biomedical database server if required.
    :param str date_field: projects deposited in PRIDE are search based on date, either
        submissionData or publicationDate (default)
    """
    ftp_pride = (  # same as data ["_links"]["datasetFtpUrl"]["href"]
        "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/YEAR/MONTH/PXDID/FILE_NAME"
    )
    url_pride_api = "http://www.ebi.ac.uk/pride/ws/archive/v3/projects/" + pxd_id
    data = None
    r = requests.get(url_pride_api, timeout=60)
    data = r.json()
    submission_date = data[date_field]
    year, month, _ = submission_date.split("-")

    ftp_url = (
        ftp_pride.replace("YEAR", year)
        .replace("MONTH", month)
        .replace("PXDID", pxd_id)
        .replace("FILE_NAME", file_name)
    )
    download_from_ftp(ftp_url, user, password, to, file_name)

    return data
