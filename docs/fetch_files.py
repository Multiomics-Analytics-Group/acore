import re
from pathlib import Path, PosixPath
from typing import Optional

import requests

OUTPUT_URL = (
    "https://github.com/bigbio/quantms/raw/refs/heads/README_links/docs/output.md"
)


def download_file(url: str, save_path: Path, timeout: int = 20):
    """
    Download a file from a URL and save it to the specified path.

    Args:
        url (str): URL of the file to download
        save_path (Path): Path where the downloaded file will be saved
        timeout (int): Timeout for the download request in seconds

    Returns:
        bool: True if download was successful, False otherwise
    """

    # Send GET request
    response = requests.get(url, stream=True, timeout=timeout)
    response.raise_for_status()  # Raise exception for HTTP errors

    save_path = Path(save_path)
    # Create directories if they don't exist
    save_path.parent.mkdir(parents=True, exist_ok=True)

    # Write content to file
    with open(save_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)

    print(f"Successfully downloaded: {url} to {save_path}")
    return True


def replace_local_links_with_github_links(content: str, repo_url_base: str):
    """
    Replace local links in the content with links to the GitHub repository.
    The regex pattern has two parts
    1. \[(.*?)\] matches the link text.
    2. \((?!https?://)([^)]+)\) matches the relative path
                                (not starting with http/https).

    Replace with main branch URL of the repository.
    """

    res = re.sub(
        r"\[(.*?)\]\((?!https?://)([^)]+)\)",
        f"[\\1]({repo_url_base}/blob/main/\\2)",
        content,
    )
    return res


def download_and_patch_output(
    repo_url_base: str,
    file_path_in_repo: str,
    output_path: Optional[str] = None,
    insert_origin_line: bool = True,
    insert_origin_line_at: int = 1,
):
    """Download the output file from the specified URL and save it locally to use
    for Sphinx documentation purposes.

    Parameters
    ----------
    repo_url_base : str
        The base URL of the repository, using the raw format (e.g., GitHub).
    file_path_in_repo : str
        The file path within the repository. Filename is the default output name.
    output_path : Optional[str], optional
        Path where the output file will be saved, by default None
    insert_origin_line : bool, optional
        Whether to insert a line indicating the origin of the file, by default True
    insert_origin_line_at : int, optional
        The line number at which to insert the origin line, by default 1
    """
    # new folders would need to be created if necessary

    if output_path is None:
        output_path = Path(file_path_in_repo).name

    url = PosixPath(repo_url_base) / "raw/refs/heads/main" / file_path_in_repo
    url = str(url)
    url = url.replace("https:/github.com", "https://github.com")
    print(f"Downloading file from: {url}")

    assert download_file(url, output_path)  # prints success message

    with open(output_path, "r", encoding="utf-8") as f:
        content = f.readlines()

    assert len(content) > 0, "Downloaded file is empty"

    # Add hint of the origin of the file
    origin_msg = f"> file downloaded from [source]({url})"
    if insert_origin_line:
        content.insert(insert_origin_line_at, f"{origin_msg}   \n")

    content = "".join(content)  # Convert list to string

    # Replace local links with GitHub links
    content = replace_local_links_with_github_links(content, repo_url_base)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)


def download_notebook(
    repo_url_base: str,
    file_path_in_repo: str,
    output_path: Optional[str] = None,
):
    """Download a notebook file from the specified URL and save it locally.

    Parameters
    ----------
    repo_url_base : str
        The base URL of the repository, using the raw format (e.g., GitHub).
    file_path_in_repo : str
        The file path within the repository. Filename is the default output name.
    output_path : Optional[str], optional
        Path where the output file will be saved, by default None
    """
    if output_path is None:
        output_path = Path(file_path_in_repo).name

    url = PosixPath(repo_url_base) / "raw/refs/heads/main" / file_path_in_repo
    url = str(url)
    url = url.replace("https:/github.com", "https://github.com")
    print(f"Downloading notebook from: {url}")

    assert download_file(url, output_path)  # prints success message


if __name__ == "__main__":
    download_notebook(
        repo_url_base="https://github.com/biosustain/dsp_course_proteomics_intro",
        file_path_in_repo="2_data_analysis.ipynb",
        output_path="downloaded/PXD040621_analysis.ipynb",
    )
