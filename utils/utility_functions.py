from os.path import basename, splitext, join
from os import getcwd, makedirs
import gzip 
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
        makedirs(tmp_dir, exist_ok=True)

    if extension:
        full_path = join(new_directory, f"{base_name}_{extension}")
    else:
        full_path = join(new_directory, f"{base_name}")
    return full_path, tmp_dir

def decompress_gzip_if_needed(file_path):
    """
    Decompresses a gzip file if needed and returns the path to the decompressed file.

    Parameters:
        file_path (str): The path to the input file.

    Returns:
        tuple: A tuple containing:
            - str: The path to the decompressed file (or original file if not gzipped).
            - bool: True if the file was gzipped and decompressed, False otherwise.
    """
    if file_path.endswith('.gz'):
        print(f"Unzipping {file_path}")
        decompressed_path = file_path[:-3]  # Remove .gz extension
        with gzip.open(file_path, 'rb') as f_in:
            with open(decompressed_path, 'wb') as f_out:
                f_out.write(f_in.read())
        print(f"Unzipped file saved as {decompressed_path}")
        return decompressed_path, True
    else:
        return file_path, False