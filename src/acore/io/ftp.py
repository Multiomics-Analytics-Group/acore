import ftplib
from pathlib import Path


def download_from_ftp(
    ftp_url: str, user: str, password: str, to: str, file_name
) -> str:
    """Download a file from an FTP server."""
    to = Path(to)
    to.mkdir(parents=True, exist_ok=True)
    try:
        domain = ftp_url.split("/")[2]
        ftp_file = "/".join(ftp_url.split("/")[3:])
        with ftplib.FTP(domain) as ftp:
            ftp.login(user=user, passwd=password)
            with open(to / file_name, "wb") as fp:
                ftp.retrbinary("RETR " + ftp_file, fp.write)
    except ftplib.error_reply as err:
        raise ftplib.error_reply(
            "Exception raised when an unexpected reply is received from the server."
            f" {err}.\nURL:{ftp_url}"
        ) from err
    except ftplib.error_temp as err:
        raise ftplib.error_temp(
            "Exception raised when an error code signifying a temporary error."
            f" {err}.\nURL:{ftp_url}"
        ) from err
    except ftplib.error_perm as err:
        raise ftplib.error_perm(
            "Exception raised when an error code signifying a permanent error."
            f" {err}.\nURL:{ftp_url}"
        ) from err
    except ftplib.error_proto as err:
        raise ftplib.error_proto(
            "Exception raised when a reply is received from the server that does not"
            " fit the response specifications of the File Transfer Protocol."
            f" {err}.\nURL:{ftp_url}"
        ) from err
