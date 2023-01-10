import scanpy as sc

import anndata2ri
anndata2ri.activate()

import rpy2.robjects as ro
ro.numpy2ri.activate()

from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr


def deseq2(adata, formula, output_name='deseq2_output'):
    '''
    Wrapper function of DESeq2 R package, that compute differential expression by condition (formula).
    ----------
    adata : AnnData
        Annotated data matrix.
    formula : string
        Formula for DESeq2 model
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

    # Save DESeq2 output
    save_deseq2_results(adata, output_name)


def save_deseq2_results(adata, output_name):
    '''
    This saves AnnData Object with DESeqDataSet object into two file:
    .rds file (DESeqDataSet object) and
    .h5ad file (AnnData object). 
    ----------
    adata : AnnData
        Annotated data matrix.
    file_name : string
        Save file name.
    '''

    # import saveRDS R library into Python
    saveRDS = ro.r['saveRDS']

    # Save R DESeqDataSet object with DESeq2 results into .rds file
    saveRDS(adata.uns["deseq2"], f'{output_name}.rds')

    del adata.uns["deseq2"]

    # Save Python object AnnData object without DESeq2 results into .h5ad file
    adata.write(f'{output_name}.h5ad')


def read_deseq2_output(output_name):
    '''
    This reads DESeqDataSet object and AnnData object and returns AnnData Object with DESeq2 results.
    ----------
    adata : AnnData
        Annotated data matrix.
    file_name : string
        Save file name.
    '''

    # Read .h5ad file with AnnData object
    adata = sc.read(f'{output_name}.h5ad')

    # Import readRDS library
    readRDS = ro.r['readRDS']

    # Read .rds file and add it into AnnData object
    adata.uns["deseq2"] = readRDS(f'{output_name}.rds')

    return adata