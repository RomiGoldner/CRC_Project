import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from adjustText import adjust_text
from typing import Callable, Optional, Tuple


SAVEFIG_DPI = 1200

# create Anndata object to use for visualization
def create_anndata(embeddings, sequences_df, n_comps=50):
    anndata = AnnData(embeddings, obs=sequences_df)
    sc.pp.pca(anndata, n_comps=n_comps)
    sc.pp.neighbors(anndata)
    sc.tl.umap(anndata)
    return anndata


# plot embedding in UMAP
def plot_embedding(color_embed, color_map, title=None, legend_size=7, n_comps=50, plot_pdf_path=None, anndata=None, fname=None):
    anndata = create_anndata(n_comps) if anndata is None else anndata
    fig = plot_anndata_rep(
        anndata,
        color=color_embed,
        direct_label=False,
        cmap=color_map,
        title=title,
        legend_size=legend_size,
        fname=fname,
    )
    # use title connected with underscore as file name
    fname = fname if fname is not None else title.replace(" ", "_") + ".svg"
    print("Saving figure to {}".format(fname))
    if fname is not None:
        fig.savefig(fname, bbox_inches="tight", dpi=SAVEFIG_DPI)
    if fig is not None:
        plt.show()
    

# plot anndata representation
def plot_anndata_rep(
    a: AnnData,
    color: str,
    representation: str = "umap",
    representation_axes_label: str = "",
    swap_axes: bool = False,
    sort_column: Optional[str] = None,
    cmap: Callable = plt.get_cmap("tab20"),
    direct_label: bool = True,
    adjust: bool = False,
    ax_tick: bool = False,
    title: str = "Title of Figure",
    legend_size: Optional[int] = None,
    figsize: Tuple[float, float] = (6.4, 4.8),
    fname: str = "",
    **kwargs,
):
    """
    Plot the given adata's representation, directly labelling instead of using
    a legend
    """
    def sort_priority(label):
        return {
            'MN_MAIT': 2,
            'Other_MAIT': 1,
            'Not_MAIT': 0
        }.get(label, -1)  # Default for any unexpected category

    if sort_column and sort_column in a.obs.columns:
        # Sort 'a.obs' according to the custom priority
        a.obs['sort_priority'] = a.obs[sort_column].apply(sort_priority)
        sorted_indices = np.argsort(a.obs['sort_priority'].values)
        a = a[sorted_indices]

        
    rep_key = "X_" + representation
    assert (
        rep_key in a.obsm
    ), f"Representation {representation} not fount in keys {a.obsm.keys()}"

    coords = a.obsm[rep_key]
    if swap_axes:
        coords = coords[:, ::-1]  # Reverse the columns
    assert isinstance(coords, np.ndarray) and len(coords.shape) == 2
    assert coords.shape[0] == a.n_obs
    assert color in a.obs
    color_vals = a.obs[color]
    unique_val = np.unique(color_vals.values)
    #unique_val = color_vals.values.categories.values # for appearance bins
    color_idx = [sorted(list(unique_val)).index(i) for i in color_vals]
    #color_idx = [list(unique_val).index(i) for i in color_vals] # for appearance bins
    # Vector of colors for each point
    color_vec = [cmap(i) for i in color_idx]

    fig, ax = plt.subplots(figsize=figsize, dpi=300)
    ax.scatter(
        coords[:, 0], coords[:, 1], s=12000 / coords.shape[0], c=color_vec, alpha=0.9
    )

    if direct_label:
        # Label each cluster
        texts = []
        for v in unique_val:
            v_idx = np.where(color_vals.values == v)
            # Median behaves better with outliers than mean
            v_coords = np.median(coords[v_idx], axis=0)
            t = ax.text(
                *v_coords,
                v,
                horizontalalignment="center",
                verticalalignment="center",
                size=legend_size,
            )
            texts.append(t)
        if adjust:
            adjust_text(
                texts,
                only_move={"texts": "y"},
                force_text=0.01,
                autoalign="y",
                avoid_points=False,
            )
    else:
        patches = []
        for i, val in enumerate(unique_val):
            p = mpatches.Patch(color=cmap(i), label=val)
            patches.append(p)
        ax.legend(handles=patches, prop={"size": legend_size})

    rep_str = representation_axes_label if representation_axes_label else representation
    if not swap_axes:
        ax.set(
            xlabel=f"{rep_str.upper()}1", ylabel=f"{rep_str.upper()}2")
    else:
        ax.set(
            xlabel=f"{rep_str.upper()}2", ylabel=f"{rep_str.upper()}1")
    ax.title.set_text(title)
    ax.set(**kwargs)
    if not ax_tick:
        ax.set(xticks=[], yticks=[])

    if fname:
        fig.savefig(fname, bbox_inches="tight", dpi=SAVEFIG_DPI)
    
    return fig


# Generate combined palette of the given number of colors 
def generate_combined_palette(n_colors):
    # Start with tab20 series which gives 60 colors, converting tuples to a list
    colors = list(plt.get_cmap('tab20').colors) + list(plt.get_cmap('tab20b').colors) + list(plt.get_cmap('tab20c').colors)
    # Need additional colors, sample from 'viridis'
    additional_needed = n_colors - len(colors)
    if additional_needed > 0:
        additional_colors = plt.get_cmap('viridis')
        # Ensure even spacing of additional colors across the colormap
        additional_indices = [i / additional_needed for i in range(additional_needed)]
        colors += [additional_colors(i) for i in additional_indices]
    return colors[:n_colors]


# plot anndata representation with direct label
def plot_anndata_rep_direct_label(
    a: AnnData,
    color: str,
    label: str,
    representation: str = "umap",
    representation_axes_label: str = "",
    swap_axes: bool = False,
    sort_column: Optional[str] = None,
    cmap: Callable = plt.get_cmap("tab20"),
    direct_label: bool = True,
    adjust: bool = False,
    ax_tick: bool = False,
    title: str = "Title of Figure",
    legend_size: Optional[int] = None,
    figsize: Tuple[float, float] = (6.4, 4.8),
    fname: str = "",
    **kwargs,
):
    """
    Plot the given adata's representation, directly labeling instead of using
    a legend.
    """
    if sort_column and sort_column in a.obs:
        # Sort the AnnData object by the specified column
        a = a[np.argsort(a.obs[sort_column])]
        
    rep_key = "X_" + representation
    assert (
        rep_key in a.obsm
    ), f"Representation {representation} not found in keys {a.obsm.keys()}"

    coords = a.obsm[rep_key]
    if swap_axes:
        coords = coords[:, ::-1]  # Reverse the columns
    assert isinstance(coords, np.ndarray) and len(coords.shape) == 2
    assert coords.shape[0] == a.n_obs
    assert color in a.obs
    assert label in a.obs
    color_vals = a.obs[color]
    label_vals = a.obs[label]
    unique_val = np.unique(color_vals.values)
    color_idx = [sorted(list(unique_val)).index(i) for i in color_vals]
    color_vec = [cmap(i) for i in color_idx]

    fig, ax = plt.subplots(figsize=figsize, dpi=300)
    ax.scatter(
        coords[:, 0], coords[:, 1], s=12000 / coords.shape[0], c=color_vec, alpha=0.9
    )

    if direct_label:
        # Label each point with the Subject
        texts = []
        for i, (x, y) in enumerate(coords):
            t = ax.text(
                x, y,
                label_vals.iloc[i],
                horizontalalignment="center",
                verticalalignment="center",
                size=legend_size - 2,
            )
            texts.append(t)
        if adjust:
            adjust_text(
                texts,
                only_move={"texts": "y"},
                force_text=0.01,
                autoalign="y",
                avoid_points=False,
            )
    else:
        patches = []
        for i, val in enumerate(unique_val):
            p = mpatches.Patch(color=cmap(i), label=val)
            patches.append(p)
        ax.legend(handles=patches, prop={"size": legend_size})

    rep_str = representation_axes_label if representation_axes_label else representation
    if not swap_axes:
        ax.set(
            xlabel=f"{rep_str.upper()}1", ylabel=f"{rep_str.upper()}2")
    else:
        ax.set(
            xlabel=f"{rep_str.upper()}2", ylabel=f"{rep_str.upper()}1")
    ax.title.set_text(title)
    ax.set(**kwargs)
    if not ax_tick:
        ax.set(xticks=[], yticks=[])

    if fname:
        fig.savefig(fname, bbox_inches="tight", dpi=SAVEFIG_DPI)
    
    return fig
