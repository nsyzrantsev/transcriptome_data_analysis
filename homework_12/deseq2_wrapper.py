import scanpy as sc

import anndata2ri
anndata2ri.activate()

import rpy2.robjects as ro
ro.numpy2ri.activate()

from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr


def deseq2(adata, formula):
    '''
    Wrapper function of DESeq2 R package, 
    that compute differential expression 
    with single cell data.
    ----------
    adata : AnnData
        Annotated data matrix.
    formula : str
        Formula for DESeq2 model

    Returns
    ----------
    AnnData object.
    '''

    # Import DESeq2 library from R to Python
    deseq2 = importr("DESeq2")

    # Create DESeqDataSet from AnnData object
    dds = deseq2.DESeqDataSetFromMatrix(countData=adata.X.T, 
                                        colData=adata.obs,
                                        design=Formula(formula))
    
    # Evaluate differential expression by DESeq2 model
    dds = deseq2.DESeq(dds)

    # Write DESeq2 results into AnnData object    
    adata.uns["deseq2"] = dds

    return adata


def write_deseq2(adata, output):
    '''
    This saves AnnData Object with DESeqDataSet object into two file:
    .rds file (DESeqDataSet object) and
    .h5ad file (AnnData object). 
    ----------
    adata : AnnData
        Annotated data matrix.
    output : str
        Save file name.
    '''

    # import saveRDS R library into Python
    saveRDS = ro.r['saveRDS']

    # Save R DESeqDataSet object with DESeq2 results into .rds file
    saveRDS(adata.uns["deseq2"], f'{output}.rds')

    del adata.uns["deseq2"]

    # Save Python object AnnData object without DESeq2 results into .h5ad file
    adata.write(f'{output}.h5ad')


def read_deseq2(read_name):
    '''
    This reads DESeqDataSet object and AnnData object
    and returns AnnData Object with DESeq2 results.
    ----------
    adata : AnnData
        Annotated data matrix.
    read_name : str
        Read file name.

    Returns
    ----------
    AnnData object.
    '''

    # Read .h5ad file with AnnData object
    adata = sc.read_h5ad(f'{read_name}.h5ad')

    # Import readRDS library
    readRDS = ro.r['readRDS']

    # Read .rds file and add it into AnnData object
    adata.uns["deseq2"] = readRDS(f'{read_name}.rds')

    return adata


def deseq2_clusters(adata, 
                    formula, 
                    clusters_column='leiden', 
                    output='deseq2_clusters'
                    ):
    '''
    Wrapper function of DESeq2 R package, 
    that compute differential expression 
    with clusterized single cell data.
    ----------
    adata : AnnData
        Annotated data matrix.
    formula : str
        Formula for DESeq2 model
    clusters_column : str
        Clusters column name.
    output : str
        Save file name.
    
    Returns
    ----------
    AnnData object.
    '''

    # Create list of cluster indexes
    clusters = adata.obs[clusters_column].value_counts().index.astype(str)
    
    # Intialized dictionary for saving 
    # DESEq2 results for every cluster
    adata.uns["deseq2"] = dict()

    # Evaluate differential expressions for every cluster
    for cluster in clusters:

        # Select AnnData only for one cluster
        cluster_adata = adata[adata.obs[clusters_column] == cluster]

        # Run DESeq2 for current cluster
        adata_with_dds = deseq2(cluster_adata, formula)

        # Extract DESeq2 results
        dds = adata_with_dds.uns["deseq2"]

        # Save DESeq2 results into dictionary
        adata.uns["deseq2"][cluster] = dds

    return adata