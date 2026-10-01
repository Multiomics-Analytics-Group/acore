import rarfile


# ! unrar must be installed in the system
def unrar(filepath, to):
    """
    Decompress RAR file
    :param str filepath: path to rar file
    :param str to: where to extract all files

    """
    try:
        with rarfile.RarFile(filepath) as opened_rar:
            print(f"Extracting files: {opened_rar.namelist()}")
            opened_rar.extractall(to)
    except Exception as err:
        print(f"Error: {filepath}. Could not unrar file {err}")
        raise
