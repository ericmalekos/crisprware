from os.path import basename, splitext, join
from os import getcwd, makedirs

def create_output(file_path, outdir=None, extension="", stripped="", tmp=False):
    """
    Creates an output directory based on the given output prefix and base directory.
    
    Parameters:
    - file_path: Path of the input file.
    - extension: Suffix to append to the directory and file name.
    - stripped: String to remove from the base_name before proceeding.
    - tmp: Boolean to create a temporary subdirectory.

    Returns:
    - full_path: The final path to the created output file.
    - tmp_dir: Path to the temporary directory (if tmp=True), otherwise an empty string.
    """

    # Extract the base name without extension
    base_name = splitext(basename(file_path))[0]

    if stripped:
        base_name = base_name.replace(stripped, "")

    # Get the current working directory
    print(outdir)
    if not outdir:
        outdir = getcwd()

    # Create a new directory name with the suffix in the current directory
    if extension:
        new_directory = join(outdir, f"{base_name}_{extension}")
    else:
        new_directory = join(outdir, f"{base_name}")

    # Create the directory
    makedirs(new_directory, exist_ok=True)

    tmp_dir = ""
    if tmp:
        tmp_dir = join(new_directory, "tmp/")
        print(new_directory)
        print(tmp_dir)
        makedirs(tmp_dir, exist_ok=True)

    if extension:
        full_path = join(new_directory, f"{base_name}_{extension}")
    else:
        full_path = join(new_directory, f"{base_name}")
    return full_path, tmp_dir