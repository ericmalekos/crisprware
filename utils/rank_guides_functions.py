import pandas as pd
from pybedtools import BedTool

from utils.gtf_bed_processing_functions import extract_ids

def create_combined_weighted_column(df, column_names, weights=None, normalize_columns = True):
    # Check if weights are provided, else set to 1 for each column
    if weights is None or not weights:
        weights = [1] * len(column_names)

    # Ensure the lists of columns and weights match in length
    if len(column_names) != len(weights):
        raise ValueError("The length of column_names and weights must match.")


    # Apply Min-Max normalization to specified columns and add as new columns
    new_col_names = []
    if normalize_columns:
        for col in column_names:
            new_col_name = f'{col}_normalized'
            new_col_names.append(new_col_name)
            df[new_col_name] = ((df[col] - df[col].min()) / (df[col].max() - df[col].min())).round(4)


    # Calculate the weighted sum
    df['combined_weighted'] = (sum(df[col] * weight for col, weight in zip(new_col_names, weights)) / len(new_col_names)).round(4)

    return df

def filter_df_by_column(df, column, min_value):
    print(f"\tFiltering by {column}, cut-off:   {min_value}")
    df = df[df[column] >= min_value]
    return df

def group_and_minimize(df, rank_column, num_to_keep):
    df = df.sort_values(by=['target_id', rank_column], ascending=[True, False])
    if num_to_keep == -1:
        print("\t--num_to_keep=-1, all gRNAs returned")
        return df

    print(f"\tKeeping top {num_to_keep} gRNAs per target, according to {rank_column}")
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

# faster implementation for spacing guides?
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
    df.columns = header

    #if the fourth column entries are not unique
    if not df.iloc[:, 3].nunique() == len(df):\
        df['target_id'] = df.iloc[:, 0].astype(str) + "_" + df.iloc[:, 1].astype(str) + "_" + df.iloc[:, 2].astype(str)
    df.columns = header

    target_ids = set(df["target_id"])

    return df, target_ids


def gRNA_to_bed(gRNAs, targets, header, target_header):

    df = pd.read_table(BedTool(gRNAs).intersect(targets, wo=True).fn, sep="\t", header=None)

    #print(df)
    #drop the last column (number of overlapping bases)
    df = df.drop(df.columns[-1], axis=1)

    num_columns = len(df.columns)

    gRNAheader = header + target_header
    columns_to_drop = target_header[:-1]

    # Determine the number of columns to rename (minimum of the two lengths)
    num_rename = min(num_columns, len(gRNAheader))

    # Create a dictionary mapping current column names to new names
    rename_dict = dict(zip(df.columns[:num_rename], gRNAheader[:num_rename]))

    # Rename the columns
    df.rename(columns=rename_dict, inplace=True)
    df.drop(columns=columns_to_drop, inplace=True)

    target_ids = set(df["target_id"])

    return df, target_ids

def gRNA_to_tscript(gRNAs, mode, targets, header):

    df = pd.read_table(BedTool(gRNAs).intersect(targets, wo=True).fn, sep="\t", header=None)
    df.columns = header + df.columns[len(header):].tolist()
    # This is a hacky way to determine t
    int_columns = [col for col in df.columns if isinstance(col, int)]
    target_col = int_columns.pop(-2)
    df = df.drop(int_columns, axis=1)
    #print(df.head())

    transcript_pattern = r'transcript_id[= ]"?([^";]*)"?'
    gene_pattern = r'gene_id[= ]"?([^";]*)"?'

    df['transcript_id'] = df[target_col].str.extract(transcript_pattern)
    df['gene_id'] = df[target_col].str.extract(gene_pattern)

    if mode == "tx":
        df.rename(columns={'transcript_id': 'target_id'}, inplace=True)
    else:
        df = df.drop(['transcript_id',target_col], axis=1)
        df.rename(columns={'gene_id': 'target_id'}, inplace=True)

    target_ids = set(df["target_id"])
    return df, target_ids

def analyze_target_ids(df, no_gRNASet):
    """
    Analyzes the distribution of gRNAs across target IDs in a DataFrame and provides statistics.

    This function calculates and prints the median, minimum, and maximum number of gRNAs associated
    with target IDs in the provided DataFrame. It also reports the number of targets with no associated
    gRNA guides based on an external set of target IDs.

    Parameters:
    - df (pd.DataFrame): DataFrame containing gRNA data with a 'target_id' column.
    - no_gRNASet (set): A set of target IDs that have no associated gRNA guides.

    Returns:
    - pd.Series: A Series containing the counts of gRNAs per target ID.

    Prints:
    - Median number of gRNAs per target.
    - Minimum number of gRNAs per target.
    - Maximum number of gRNAs per target.
    - Number of targets with 0 gRNA guides.
    """

    # Count the occurrences of each target_id and calculate basic statistics
    counts = df['target_id'].value_counts()

    zero_counts = pd.Series(0, index=no_gRNASet)
    all_counts = pd.concat([counts, zero_counts])
    median_all = all_counts.median()


    median_count = counts.median()
    min_count = counts.min()
    max_count = counts.max()

    print(f"\n#\tMedian number of gRNAs per target: {median_all}")
    print(f"#\tNumber of targets with 0 gRNA guides: {len(no_gRNASet)}\n#")


    print(f"#\tCalculations exluding targets with 0 counts:")
    print(f"#\n#\t\tMedian number of gRNAs per target: {median_count}")
    print(f"#\t\tMinimum number of gRNAs per target: {min_count}")
    print(f"#\t\tMaximum number of gRNAs per target: {max_count}")
    print("#\n#######################################################################\n")
    # Return the counts of gRNAs per target ID for further analysis if needed
    return counts
