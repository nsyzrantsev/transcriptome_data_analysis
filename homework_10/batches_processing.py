import sys

# Libraries for data analysis
import scanpy as sc
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

# Load batch correction method
import harmonypy

# Load libraries for working with R
import rpy2.robjects as ro
import anndata2ri

ro.numpy2ri.activate()
anndata2ri.activate()

def droplet_clean(adata):
    '''
    This clean AnnData object from droplets and select statistical significant cells.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.

    Returns
    ----------
    AnnData object.
    '''
    adata = adata[adata.X.sum(axis=1).A.T[0] >= 10]
    # load DropletUtils package in R
    ro.r("suppressPackageStartupMessages(library(DropletUtils))")
    # create input_h5ad variable in R
    ro.globalenv["input_h5ad"] = adata.copy()
    # clearing the dataset from droplets by emptyDrops function
    adata_without_droplets = ro.r('cleared_expressions <- emptyDrops(assay(input_h5ad))')
    # Select statistically significant cells
    cell_barcodes = adata_without_droplets.loc[adata_without_droplets.FDR < 0.05]
    adata_cleared = adata[cell_barcodes.index]
    return adata_cleared

def select_hvg(adata, batch_key):
    '''
    This selects high-variable genes in the AnnData object by function in the scanpy package.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.

    Returns
    ----------
    AnnData object.
    '''
    adata_hvg = adata.copy()

    sc.pp.normalize_total(adata_hvg)
    sc.pp.log1p(adata_hvg)

    sc.pp.highly_variable_genes(
        adata_hvg,
        n_top_genes=3000,
        flavor="seurat_v3",
        layer="counts",
        batch_key=batch_key
    )
    adata_hvg.raw = adata_hvg
    adata_hvg = adata_hvg[:, adata_hvg.var.highly_variable]

    return adata_hvg

def add_batch_name_column(adata, batch_key, batch_value):
    '''
    This adds a batch column with values into the AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        Name of column in `adata.obs` that corresponds to the name of samples (batches).
    batch_value : str
        Cell value in `adata.obs` allows for determining whether a cell belongs to some sample (batch).

    Returns
    ----------
    AnnData object.
    '''
    adata.obs[batch_key] = [batch_value] * len(adata.obs)
    return adata

def split_batches(adata, batch_key):
    '''
    This splits one AnnData object into few ones (per batch).
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        Name of column in `adata.obs` that corresponds to the name of samples / batches.

    Returns
    ----------
    List of AnnData objects.
    '''
    adatas = []
    for batch in set(adata.obs[batch_key]):
        adatas.append(adata[adata.obs[batch_key] == batch].copy())
    return adatas

def concatenate_batches(adatas, batch_key="batch"):
    '''
    This function concatenates different AnnData objects into one.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str, optional
        Name of column in `adata.obs` that corresponds to the name of samples / batches.

    Returns
    ----------
    AnnData object.
    '''
    adata = adatas[0].concatenate(adatas[1:], batch_key=batch_key, index_unique=None)
    return adata

def harmony(adata, batch_key):
    '''
    This evaluates harmony coordinates for visualization clusters in AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str, optional
        Name of column in `adata.obs` that corresponds to the name of samples / batches.

    Returns
    ----------
    AnnData object.
    '''
    adata_harmony = adata.copy()

    sc.external.pp.harmony_integrate(adata_harmony, key=batch_key, max_iter_harmony=20)

    sc.pp.neighbors(
        adata_harmony,
        use_rep="X_pca_harmony",
        n_pcs=30
    )

    sc.tl.leiden(adata_harmony)

    sc.tl.umap(adata_harmony)

    return adata_harmony

def filter_adata(h5ad_file, batch_key, batch_value):
    '''
    This filters the AnnData object and prepares him for the next steps.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str, optional
        Name of column in `adata.obs` that corresponds to the name of samples / batches.
    batch_value : str
        Cell value in `adata.obs` allows for determining whether a cell belongs to some sample (batch).

    Returns
    ----------
    AnnData object.
    '''
    adata = sc.read_h5ad(h5ad_file)
    adata.layers["counts"] = adata.X.copy()
    batched_adata = add_batch_name_column(adata, batch_key, batch_value)
    adata_without_droplets = droplet_clean(batched_adata)
    adata_hvg = select_hvg(adata_without_droplets, batch_key)
    return adata_hvg

if __name__ == '__main__':
    # Load batch correction method type 
    # into this script from command line argument
    h5ad_list = sys.argv[1:]
    batch_key = 'batch'

    # Filter all datasets and
    # save them into list
    adata_list = []
    for i, h5ad_file in enumerate(h5ad_list):
        batch_value = f'{batch_key}_{i+1}'
        filtered_adata = filter_adata(h5ad_file, batch_key, batch_value)
        adata_list.append(filtered_adata)

    # Concatenate adatas from different samples
    # into one concatenated adata object
    concatenated_adatas = concatenate_batches(adata_list, batch_key)

    # Scale and calculate PCA 
    # for all adatas in an adata list
    sc.pp.scale(concatenated_adatas)
    sc.tl.pca(concatenated_adatas)

    # Run harmony batch correction method
    # for concatenated batches
    batches_after_harmony_correction = harmony(concatenated_adatas, batch_key)

    # Save integrated batches as a .h5ad file
    batches_after_harmony_correction.write_h5ad("integrated_batches.h5ad")