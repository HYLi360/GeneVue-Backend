def file_renamer(filename_list: list[str], prefix: str = "", suffix: str = ""):
    return ["".join([prefix, filename, suffix]) for filename in filename_list]
