import os
import shutil
import subprocess
from pathlib import Path

import rarfile

_UNRAR_FALLBACK_DIRS = ("/opt/homebrew/bin", "/usr/local/bin", "/usr/bin")


def _find_unrar() -> str | None:
    """Find the unrar binary, also in common install locations that may be missing
    from PATH (e.g. in a notebook kernel started from an IDE)."""
    return shutil.which("unrar") or shutil.which(
        "unrar", path=os.pathsep.join(_UNRAR_FALLBACK_DIRS)
    )


# ! unrar must be installed in the system
def unrar(filepath, to):
    """
    Decompress RAR file

    Uses the system ``unrar x`` command, as extraction through ``rarfile`` can fail
    for some archives (``rarfile`` falls back to other tools such as ``bsdtar`` if
    ``unrar`` is not on PATH). Falls back to ``rarfile`` if ``unrar`` is not found.

    :param str filepath: path to rar file
    :param str to: where to extract all files

    """
    filepath = Path(filepath).resolve()
    to = Path(to)
    to.mkdir(parents=True, exist_ok=True)
    try:
        unrar_bin = _find_unrar()
        if unrar_bin is not None:
            # -o+ overwrite existing files, -y assume yes on all queries
            subprocess.run(
                [unrar_bin, "x", "-o+", "-y", str(filepath), f"{to.resolve()}/"],
                check=True,
                capture_output=True,
                text=True,
            )
        else:
            print("`unrar` not found on PATH, falling back to rarfile.")
            with rarfile.RarFile(filepath) as opened_rar:
                print(f"Extracting files: {opened_rar.namelist()}")
                opened_rar.extractall(to)
    except subprocess.CalledProcessError as err:
        print(f"Error: {filepath}. Could not unrar file {err.stderr or err.stdout}")
        raise
    except Exception as err:
        print(f"Error: {filepath}. Could not unrar file {err}")
        raise
