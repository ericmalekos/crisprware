import pandas as pd
import os

def calculate_statistics(df):
    """
    Calculates statistical metrics for transcript expression levels in a DataFrame.
    Parameters:
    - df (pandas.DataFrame): A DataFrame containing transcript expression levels. Column names should include
                             'TPM_' or 'COUNTS_' prefixes to indicate expression measurements.

    Returns:
    - pandas.DataFrame: The original DataFrame augmented with four new columns containing the calculated
                        statistics for each transcript. The numerical results are rounded to three decimal places.
    """
    
    tpm_columns = [col for col in df.columns if 'TPM_' in col or 'COUNTS_' in col]
    
    df['tscript_min'] = df[tpm_columns].min(axis=1)
    df['tscript_max'] = df[tpm_columns].max(axis=1)
    df['tscript_median'] = df[tpm_columns].median(axis=1)
    df['tscript_mean'] = df[tpm_columns].mean(axis=1)
    return df.round(3)

def process_dataframes(dfs):
    """
    Merges multiple DataFrames and calculates statistical metrics on the merged data.

    Parameters:
    - dfs (list of pandas.DataFrame): A list of DataFrames to be merged. Each DataFrame should have the same
                                      first column which will be used as the key for merging.

    Returns:
    - pandas.DataFrame: A DataFrame containing the merged data from all input DataFrames, with statistical metrics
                        calculated across the columns. The specific metrics calculated are determined by the
                        `calculate_statistics` function.
    """
    
    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on=df.columns[0], how='outer')
    return calculate_statistics(merged_df)


def process_kallisto(files):
    """
    Processes Kallisto quantification files and aggregates or analyzes their data.
    
    Parameters:
    - files (list of str): A list of file paths to Kallisto quantification output files. Each file should be
                           a tab-separated values (TSV) file with at least 'target_id' and 'tpm' columns.

    Returns:
    - pandas.DataFrame: Depending on the number of files provided, returns either:
                        - A single DataFrame with calculated statistics for one file's data, or
                        - An aggregated DataFrame with data from multiple files, ready for further processing.
    """
    
    dfs = []
    for file in files:
        data = pd.read_csv(file, sep='\t')
        # Renaming 'target_id' column to 'transcript_id'
        dfs.append(data[['target_id', 'tpm']].rename(columns={'target_id': 'transcript_id', 'tpm': f'TPM_{file}'}))

    if len(dfs) > 1:
        return process_dataframes(dfs)
    else:
        return calculate_statistics(dfs[0])

def process_salmon(files):
    """
    Similar to process_kallisto(files)
    """
    
    dfs = []
    for file in files:
        data = pd.read_csv(file, sep='\t')
        # Renaming 'Name' column to 'transcript_id'
        dfs.append(data[['Name', 'TPM']].rename(columns={'Name': 'transcript_id', 'TPM': f'TPM_{file}'}))
    
    if len(dfs) > 1:
        return process_dataframes(dfs)
    else:
        return calculate_statistics(dfs[0])

def process_flair(files):
    """
    Similar to process_kallisto(files)
    """
    
    dfs = []
    for file in files:
        data = pd.read_csv(file, sep='\t')
        # Renaming 'target_id' column to 'transcript_id'
        data.columns = [data.columns[0]] + [f'COUNTS_{col}' for col in data.columns[1:]]
        data.iloc[:, 0] = data.iloc[:, 0].str.split('_').str[0]
        data['ids'] = data['ids'].str.replace(';', ':')
        dfs.append(data.rename(columns={'ids': 'transcript_id'}))
        print(data.head())

    if len(dfs) > 1:
        return process_dataframes(dfs)
    else:
        return calculate_statistics(dfs[0])

def process_mandalorian(files):
    """
    Similar to process_kallisto(files)
    """
    
    dfs = []
    for file in files:
        data = pd.read_csv(file, sep='\t')
        # Renaming 'target_id' column to 'transcript_id'
        data.columns = [data.columns[0]] + [f'TPM_{col}' for col in data.columns[1:]]
        dfs.append(data.rename(columns={'Isoform': 'transcript_id'}))

    if len(dfs) > 1:
        return process_dataframes(dfs)
    else:
        return calculate_statistics(dfs[0])

def infer_file_type_from_first_line(path):
    """
    Infers the file type of a given file based on its first line.

    Parameters:
    - path (str): The file path of the file to be checked.

    Returns:
    - str or None: Returns a string identifier of the file type ('kallisto', 'salmon', or 'flair') if a match is found.
                   If the file type cannot be determined or the file cannot be read, returns None.
    """

    print("\n\tInferring file type from header line\n")
    try:
        with open(path, 'r') as file:
            first_line = file.readline().strip()
        
        if first_line.startswith('target_id\tlength\teff_length\test_counts\ttpm'):
            print("\t\t" + path + " is a Kallisto file")
            return 'kallisto'
        
        elif first_line.startswith('Name\tLength\tEffectiveLength\tTPM\tNumReads'):
            print("\t\t" + path + " is a Salmon file")
            return 'salmon'
        
        elif first_line.startswith('ids\t'):
            print("\t\t" + path + " is a FLAIR file")
            return 'flair'
        
        elif first_line.startswith('Isoform\t'):
            print("\t\t" + path + " is a Mandalorian file")
            return 'mandalorian'
        
        
        else:
            print("\n\t\tCould not determine file type for " + path)
            print("\t\tCheck file header")
            return None
    
    except Exception as e:
        print(f"\n\tCould not read file {path} to infer file type: {e}")
        return None

def process_files(paths, filetype="infer"):
    """
    Processes files based on their type, handling different genomic data formats.

    Parameters:
    - paths (list of str or str): A single file path or a list of file paths to be processed.
    - filetype (str, optional): The type of the files to process ('kallisto', 'salmon', 'flair', 'mandalorian', or 'infer').
                                Default is 'infer', which means the function will try to determine the file type
                                automatically based on file content.

    Returns:
    - The result of the file processing, which depends on the specific processing function called for the detected
      file type. This could vary from aggregated data frames to statistical analysis results.

    Raises:
    - ValueError: If any of the paths is not a file, if an invalid file type is specified or inferred,
                  if multiple different file types are found in the input list, or if no valid files are found.


    """
        
    if not isinstance(paths, list):
        paths = [paths]

    kallisto_files = []
    salmon_files = []
    flair_files = []
    mandalorian_files = []

    for path in paths:
        if os.path.isfile(path):
            if filetype == "infer":
                filetype = infer_file_type_from_first_line(path)
            
            #print("\n\tProcessing files as " + filetype)
            if filetype == "kallisto":
                kallisto_files.append(path)
            elif filetype == "salmon":
                salmon_files.append(path)
            elif filetype == "flair":
                flair_files.append(path)
            elif filetype == "mandalorian":
                mandalorian_files.append(path)
            else:
                raise ValueError(f"Invalid file type for path: {path}")
        else:
            raise ValueError(f"Path is not a file: {path}")

    if sum([bool(kallisto_files), bool(salmon_files), bool(flair_files)]) == 0:
        raise ValueError("Unable to infer filetype")
    if sum([bool(kallisto_files), bool(salmon_files), bool(flair_files)]) > 1:
        raise ValueError("Multiple filetypes passed. Please pass only one of Salmon/Kallisto/FLAIR")
  
    if kallisto_files:
        return process_kallisto(kallisto_files)
    elif salmon_files:
        return process_salmon(salmon_files)
    elif flair_files:
        return process_flair(flair_files)
    elif mandalorian_files:
        return process_mandalorian(mandalorian_files)
    else:
        raise ValueError("No valid files found.")

def filter_dataframe(df, **filters):
    """
    Filters a DataFrame based on specified column thresholds.

    Parameters:
    - df (pandas.DataFrame): The DataFrame to be filtered. It should contain the columns specified in the filters.
    - **filters (dict): Keyword arguments where keys are column names in the DataFrame and values are the minimum
                        threshold values for those columns. Rows in the DataFrame where the column value is below
                        the threshold will be excluded from the result.

    Returns:
    - pandas.DataFrame: A DataFrame containing only the rows that meet all the specified filtering criteria.

    Prints:
    - The function prints the number of unique transcripts (or entries based on the 'transcript_id' column) before 
      and after filtering to provide an overview of how the filtering criteria impacted the dataset size.
    """

    print(f'\n\tInitial unique transcripts:\t\t\t{    df["transcript_id"].nunique()}')

    for column, threshold in filters.items():
        if column in ["tscript_min", "tscript_max", "tscript_median", "tscript_mean"]:
            df = df[df[column] > threshold]
    
    print(f'\tTranscripts after filtering by expression:\t{df["transcript_id"].nunique()}')

    return df

def add_gene_ids_and_subset(df, relationship, col = 'median', n=-1):
    """
    Maps transcript IDs to gene IDs, filters, and subsets the DataFrame based on specified criteria.

    Parameters:
    - df (pandas.DataFrame): The DataFrame containing transcript data. Must include a 'transcript_id' column.
    - relationship (dict): A dictionary mapping transcript IDs (keys) to gene IDs (values).
    - col (str, optional): The column name to use for filtering transcripts within each gene. The default is 'median'.
                            This column name is prefixed with 'tscript_' in filter_dataframe() function.
    - n (int, optional): The number of top transcripts to retain per gene based on the `col` values. If `n` is -1 (default),
                         all transcripts are retained.

    Returns:
    - pandas.DataFrame: A DataFrame that has been optionally subsetted based on the provided criteria. If any transcripts
                        do not have a corresponding gene ID in the provided relationship dictionary, a warning is printed,
                        and those transcripts are excluded from the result.

    Prints:
    - Warnings about any transcripts not found in the provided relationship dictionary.
    - Information about the number of unique genes and transcripts before and after filtering.
    """
    
    # Map transcript_id to gene_id using the relationship dictionary
    df['gene_id'] = df['transcript_id'].map(relationship)

    df = df[['gene_id'] + [col for col in df if col != 'gene_id']]

    # check for duplicates
    # duplicated_entries = df[df['gene_id'].duplicated(keep=False)]    
    
    missing_transcripts = df[df['gene_id'].isna()]['transcript_id'].tolist()
    if missing_transcripts:
        print('\n\tWarning: Transcripts not found in GTF/GFF:\n\n\t' + '\n\t'.join(missing_transcripts[:10]) + '\n\t...\n\tand ' + \
                                                                         str(len(missing_transcripts) - 10) + ' more\n' )
        print('\tGENE LEVEL FILTERING OF TOP TRANSCRIPTS WILL BE IGNORED')
        print('\tCHECK THAT THE SAME ANNOTATION WAS USED FOR QUANTIFICATION AND CURRENT PROCESSING')
        print('\tCHECK THAT THE FIRST COLUMN OF THE TPM QUANTIFICATION FILES CONTAINS ONLY THE TRANSCRIPT ID\n')
        
        return df.drop('gene_id', axis=1)
            
    # Sort by the specified column in descending order and then group by gene_id 
    # to select the top n transcripts

    if n != -1:
        print(f'\tRetaining top {n} transcripts per gene')

    col = 'tscript_' + col

    def get_top_n_or_all(group: pd.DataFrame) -> pd.DataFrame:
        return group.nlargest(n, col)

    if n >= 0:
        df = df.groupby('gene_id', group_keys=False).apply(get_top_n_or_all)

    print(f'\tFinal unique genes:\t\t{df["gene_id"].nunique()}')
    print(f'\tFinal unique transcripts:\t{df["transcript_id"].nunique()}')

    return df