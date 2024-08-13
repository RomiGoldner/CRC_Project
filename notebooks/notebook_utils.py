import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans


# plot countplot according to given parameters
def plot_count(df, X, hue=None, order=None, title=None, xlabel=None, legend_title=None, labels=None, rotation=0, figure_name=None, palette='Set2'):
    # Set font size to 12
    plt.rcParams.update({'font.size': 12})

    # Create a countplot with Seaborn
    plt.figure(figsize=(10, 6))
    sns.countplot(data=df, x=X, hue=hue, order=order, palette=palette)
    
    if title is not None:
        plt.title(title)
    if xlabel is not None:
        plt.xlabel(xlabel)
    plt.ylabel('Count')
    if legend_title is not None:
        plt.legend(title=legend_title, labels=labels)
    plt.xticks(rotation=rotation)
    plt.tight_layout()
    if figure_name is not None:
        plt.savefig(figure_name, dpi=1200)
    plt.show()


# plot heatmap of specified columns from the given dataframe
def plot_CRC_unique_combination_heatmap(df, col1, col2, title, file_path):
    # Create a new DataFrame with only the necessary columns
    df_new = df[['Subject', col1, col2]].drop_duplicates()

    # Count the occurrences of each combination across all subjects
    aggregated_data = df_new.groupby([col1, col2]).size().reset_index(name='Count')

    # Pivot the data to create a matrix format suitable for heatmap visualization
    pivot_table = aggregated_data.pivot(col1, col2, 'Count')

    # set font size to 12
    plt.rcParams.update({'font.size': 18})
    
    # Visualize with a heatmap
    plt.figure(figsize=(10, 7))
    sns.heatmap(pivot_table, annot=True, cmap='coolwarm', fmt='g')
    plt.title(title)
    plt.ylabel(f'{col1} levels')
    plt.xlabel(f'{col2} levels')
    plt.tight_layout()
    plt.savefig(file_path, dpi=1200)
    plt.show()


# label MAIT cell according to V/J genes
def label_MAIT(df, v_gene_tra, j_gene_tra, cell_label):
    cell_label_not = 'non-'+cell_label
    df_tra = df[df['chain'] == 'TRA']
    df_trb = df[df['chain'] == 'TRB']
    # create "MAIT_cell" column if the v_gene is "TRAV1-2" and j_gene is TRAJ33 or TRAJ20 or TRAJ12
    df_tra[cell_label] = np.where((df_tra['v_gene'] == v_gene_tra) &
                                     (df_tra['j_gene'].str.contains('|'.join(j_gene_tra))), cell_label, cell_label_not)
    df_trb[cell_label] = cell_label_not
    # concatenate the two dataframes
    df_MAIT_joined = pd.concat([df_tra, df_trb])
    return df_MAIT_joined


# label single MAIT cells according to V/J genes
def label_MAIT_single_cell(sc_df, v_gene_list, j_gene_list, cell_label):
    cell_label_not = 'non-'+cell_label
    sc_df_tra = sc_df[sc_df['chain'] == 'TRA']

    # create "MAIT_cell" column if the v_gene is "TRAV1-2" and j_gene is TRAJ33 or TRAJ20 or TRAJ12
    sc_df_tra[cell_label] = np.where((sc_df_tra['v_gene'] == v_gene_list) &
                                     (sc_df_tra['j_gene'].str.contains('|'.join(j_gene_list))), cell_label, cell_label_not)

    # extract the list of barcodes from the dataframe that are MAIT cells
    true_label_barcodes = sc_df_tra[sc_df_tra[cell_label] == cell_label]['barcode'].tolist()

    # create a column in embed_data that is MAIT_cell if the barcode is in the MAIT_barcodes list, and non-MAIT_cell otherwise
    sc_df[cell_label] = np.where(sc_df['barcode'].isin(true_label_barcodes), cell_label, cell_label_not)
    return sc_df


# combine embedding to represent a patient using the specified method
def combine_embeddings(group, method='mean'):
    embeddings = group.iloc[:, :768].values
    if method == 'weighted_mean':
        proportions = group['Proportion'].values.reshape(-1, 1)
        total_proportion = proportions.sum()
        if total_proportion != 0:
            normalized_proportions = proportions / total_proportion
        else:
            normalized_proportions = proportions
        combined_embeddings = np.average(embeddings, axis=0, weights=normalized_proportions.flatten())
    elif method == 'mean':
        combined_embeddings = np.mean(embeddings, axis=0)
    elif method == 'sum':
        combined_embeddings = np.sum(embeddings, axis=0)
    elif method == 'clustering':
        kmeans = KMeans(n_clusters=4)
        kmeans.fit(embeddings)
        cluster_centroids = kmeans.cluster_centers_
        combined_embeddings = np.mean(cluster_centroids, axis=0)
    else:
        raise ValueError("Method not supported")
    return combined_embeddings



    