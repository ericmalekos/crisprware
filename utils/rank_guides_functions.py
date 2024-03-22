import pandas as pd
from pybedtools import BedTool

from utils.gtf_bed_processing_functions import extract_ids

def create_combined_weighted_column(df, column_names, weights=None):
    # Check if weights are provided, else set to 1 for each column
    if weights is None or not weights:
        weights = [1] * len(column_names)
    
    # Ensure the lists of columns and weights match in length
    if len(column_names) != len(weights):
        raise ValueError("The length of column_names and weights must match.")
    
    # Calculate the weighted sum
    df['combined_weighted'] = sum(df[col] * weight for col, weight in zip(column_names, weights)) / len(column_names)
    
    return df

def filter_df_by_column(df, column, min_value):
    print(f"\tFiltering by {column}, cut-off:   {min_value}")
    df = df[df[column] >= min_value]
    return df

def group_and_minimize(df, rank_column, num_to_keep):
    df = df.sort_values(by=['target_id', rank_column], ascending=[True, False])
    if num_to_keep == -1:
        print("\t--num_to_keep=-1, all sgRNAs returned")
        return df

    print(f"\tKeeping top {num_to_keep} sgRNAs per target, according to {rank_column}")
    # Group by target_id and keep the top entries as specified by number_of_guides
    df = df.groupby('target_id').head(num_to_keep)
    return df


def select_guides(df, rank_column, min_spacing):
    df =  df.sort_values(by=['target_id', rank_column], ascending=[True, False])

    if min_spacing == 0:
        print("\t--min_spacing=0, no filtering applied")
        return df

    print(f"\tEnforcing minimum spacing of {min_spacing} nucleotides between guides")

    def is_overlapping_or_close(new_start, new_stop, selected_guides):
        for start, stop in selected_guides:
            if (new_start <= stop + min_spacing) and (new_stop >= start - min_spacing):
                return True
        return False

    def select_from_group(group):
        selected_guides = []  # List to store tuples of (start, stop)

        for _, row in group.iterrows():
            if not is_overlapping_or_close(row['start'], row['stop'], selected_guides):
                selected_guides.append((row['start'], row['stop']))

        return group[group.apply(lambda row: (row['start'], row['stop']) in selected_guides, axis=1)]

    # Assuming df is already grouped by 'target_id'
    return df.groupby('target_id', group_keys=False).apply(select_from_group)

# might be a faster implementation for spacing guides.
# from intervaltree import Interval, IntervalTree
# def select_guides(df, min_spacing):
#     def select_from_group(group):
#         tree = IntervalTree()
#         selected_indices = []

#         for index, row in group.iterrows():
#             start, stop = row['start'], row['stop']
#             if not tree.overlaps(start - min_spacing, stop + min_spacing):
#                 tree.add(Interval(start, stop))
#                 selected_indices.append(index)

#         return group.loc[selected_indices]

#     return df.groupby('target_id', group_keys=False).apply(select_from_group)

def df_to_pybed(df):
    header = df.columns.tolist()
    bed_df = BedTool.from_dataframe(df)
    return bed_df, header

def validate_and_modify_bed(file_path):
    # Read the bed file, ignoring lines that start with "#"
    df = pd.read_csv(file_path, sep='\t', comment='#', header=None)

    # Check for minimum column requirement
    if df.shape[1] < 3:
        raise ValueError("BED must have at least 3 columns.")
    
    # Check if second and third columns are integers and validate their values
    if not (df.iloc[:, 1].dtype.kind in 'i' and df.iloc[:, 2].dtype.kind in 'i'):
        raise ValueError("Second and third columns must be integers.")
    
    if any(df.iloc[:, 1] >= df.iloc[:, 2]):
        raise ValueError("Values in the second column must be less than those in the third column.")
    
    # Drop additional columns if more than four
    if df.shape[1] > 3:
        df = df.iloc[:, :4]
    
    header = ['#target_chr', 'target_start', 'target_end', 'target_id']

    #if the fourth column entries are not unique
    if not df.iloc[:, 3].nunique() == len(df):\
        df['target_id'] = df.iloc[:, 0].astype(str) + "_" + df.iloc[:, 1].astype(str) + "_" + df.iloc[:, 2].astype(str)
    df.columns = header
    
    target_ids = set(df["target_id"])

    return df, target_ids
    

def sgRNA_to_bed(sgRNAs, targets, header, target_header):
    
    df = pd.read_table(BedTool(sgRNAs).intersect(targets, wo=True).fn, sep="\t", header=None)

    #print(df)
    #drop the last column (number of overlapping bases)
    df = df.drop(df.columns[-1], axis=1)

    num_columns = len(df.columns)

    sgRNAheader = header + target_header
    columns_to_drop = target_header[:-1]

    # Determine the number of columns to rename (minimum of the two lengths)
    num_rename = min(num_columns, len(sgRNAheader))

    # Create a dictionary mapping current column names to new names
    rename_dict = dict(zip(df.columns[:num_rename], sgRNAheader[:num_rename]))

    # Rename the columns
    df.rename(columns=rename_dict, inplace=True)
    df.drop(columns=columns_to_drop, inplace=True)

    target_ids = set(df["target_id"])

    return df, target_ids


def sgRNA_to_tscript(sgRNAs, mode, targets, header):
    
    df = pd.read_table(BedTool(sgRNAs).intersect(targets, wo=True).fn, sep="\t", header=None)
    df.columns = header + df.columns[len(header):].tolist()


    # This is a hacky way to determine the 
    int_columns = [col for col in df.columns if isinstance(col, int)]

    # add either the transcript and gene id or the gene_id as "target_id"
    if mode == "transcript":
        df[['target_id', 'gene_id']] = df.iloc[:, int_columns[8]].apply(extract_ids).apply(pd.Series)
    else:
        # assign an arbitrary number as it will be dropped later
        arb_int = 999
        df[[arb_int, 'target_id']] = df.iloc[:, int_columns[8]].apply(extract_ids).apply(pd.Series)
    # Dropping integer-labeled columns
    int_columns.append(arb_int)
    df.to_csv('test1.tsv', sep='\t', index=False, quoting=3)
    df = df.drop(int_columns, axis=1)
    df.to_csv('test2.tsv', sep='\t', index=False, quoting=3)  

    target_ids = set(df["target_id"])


    return df, target_ids

def analyze_target_ids(df, no_sgRNASet):
    """
    Analyzes the distribution of sgRNAs across target IDs in a DataFrame and provides statistics.

    This function calculates and prints the median, minimum, and maximum number of sgRNAs associated
    with target IDs in the provided DataFrame. It also reports the number of targets with no associated
    sgRNA guides based on an external set of target IDs.

    Parameters:
    - df (pd.DataFrame): DataFrame containing sgRNA data with a 'target_id' column.
    - no_sgRNASet (set): A set of target IDs that have no associated sgRNA guides.

    Returns:
    - pd.Series: A Series containing the counts of sgRNAs per target ID.

    Prints:
    - Median number of sgRNAs per target.
    - Minimum number of sgRNAs per target.
    - Maximum number of sgRNAs per target.
    - Number of targets with 0 sgRNA guides.
    """

    # Count the occurrences of each target_id and calculate basic statistics
    counts = df['target_id'].value_counts()

    median_count = counts.median()
    min_count = counts.min()
    max_count = counts.max()

    print(f"\n\tMedian number of sgRNAs per target: {median_count}")
    print(f"\tMinimum number of sgRNAs per target: {min_count}")
    print(f"\tMaximum number of sgRNAs per target: {max_count}")
    print(f"\tNumber of targets with 0 sgRNA guides: {len(no_sgRNASet)}\n")

    # Return the counts of sgRNAs per target ID for further analysis if needed
    return counts