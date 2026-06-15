"""Download files from the internet."""

import requests


def download_file(url: str, local_filename: str) -> None:
    """Download a file from the internet."""
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(local_filename, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
