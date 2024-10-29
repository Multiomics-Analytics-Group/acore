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
            print("Extracting files: {}".format(opened_rar.namelist()))
            opened_rar.extractall(to)
    except Exception as err:
        print("Error: {}. Could not unrar file {}".format(filepath, err))
        raise
