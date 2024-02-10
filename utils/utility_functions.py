from pathlib import Path


def create_output_directory(output_prefix="", base_dir="out"):
    """
    Creates an output directory based on the given output prefix and base directory.

    Parameters:
    output_prefix (str): The output prefix to append to the base directory.
    base_dir (str): The base directory for the output. Defaults to 'out/'.
    """
    full_path = base_dir + "/" + output_prefix
    full_path = full_path.replace("./", "/")
    if full_path[0] == "/": full_path = full_path[1:]
    full_path = Path(full_path.replace("//", "/"))

    if output_prefix.endswith("/"):
        full_path.mkdir(parents=True, exist_ok=True)
        full_path = str(full_path) + "/"
    else:
        full_path = str(full_path)
        ndx = full_path.rfind('/')
        dir_path=Path(full_path[:ndx])
        dir_path.mkdir(parents=True, exist_ok=True)

    return full_path