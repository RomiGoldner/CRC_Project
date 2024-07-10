import os
import pandas as pd
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def process_file_path(file_path):
    parts = file_path.split(os.sep)
    subject_type = parts[-2]  # TRA or TRB
    file_name = os.path.basename(file_path)
    subject = file_name.split('_')[0] + '_' + file_name.split('_')[1]  # e.g., pool3_S10
    return subject, subject_type


def read_and_combine_txt_files(root_dir, sep='\t'):
    """
    Reads text files from the given root directory (including subdirectories),
    extracts 'Subject' and 'Subject_Type' from each file path, and combines them into a single DataFrame.

    Parameters:
    - root_dir: The root directory to search for text files.
    - sep: The separator used in the text files (default is tab).

    Returns:
    - A pandas DataFrame containing the combined data from all text files.
    """
    all_dfs = []  # List to store individual dataframes
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith('.txt'):
                file_path = os.path.join(root, file)
                subject, subject_type = process_file_path(file_path)

                # Read the file into a DataFrame
                temp_df = pd.read_csv(file_path, sep=sep)
                temp_df['Subject'] = subject
                temp_df['Subject_Type'] = subject_type

                # Append the dataframe to the list
                all_dfs.append(temp_df)

    # Concatenate all dataframes in the list
    all_data = pd.concat(all_dfs, ignore_index=True)
    return all_data


def plot_heatmap_with_buckets(
        df, index_col='Subject', columns_col='v_gene', values_col='count', aggfunc='sum', fill_value=0, 
        percentiles=[20, 40, 60, 80, 95], figsize=(20, 10), cmap='viridis', label_col='N'):
    
    pivot_table = df.pivot_table(index=index_col, columns=columns_col, values=values_col, aggfunc=aggfunc, fill_value=fill_value)
    labels = df.groupby('Subject')[label_col].first().reindex(pivot_table.index)
    pivot_table.index = [f"(Label: {labels[idx]}): {idx} " for idx in pivot_table.index]
    pivot_table = pivot_table.sort_index()
        
    # Buckets
    pivot_values_flat = pivot_table.values.flatten()
    percentiles_values = np.percentile(pivot_values_flat, percentiles)
    bucket_thresholds = [0] + list(percentiles_values)
    def assign_bucket(value):
        for i, threshold in enumerate(bucket_thresholds):
            if value <= threshold:
                return i
        return len(bucket_thresholds)
    pivot_table_bucketed = pivot_table.applymap(assign_bucket)
    
    # Plotting
    plt.figure(figsize=figsize)
    cmap = sns.color_palette(cmap, len(bucket_thresholds) + 1)  # Create a discrete colormap
    sns.heatmap(pivot_table_bucketed, cmap=cmap, cbar_kws={'ticks': range(len(bucket_thresholds) + 1), 'label': 'Bucket'})
    plt.title(f'Heatmap of {values_col} Distribution with Percentile Buckets')
    plt.xlabel(columns_col)
    plt.ylabel(index_col)
    plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for readability
    plt.tight_layout()  # Adjust layout
    plt.show()

# Example usage:
# plot_heatmap_with_buckets(df=all_data_df_with_PP_label, index_col='Subject', columns_col='v_gene', values_col='count')
