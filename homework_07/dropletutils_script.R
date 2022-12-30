#!/usr/local/bin Rscript

library(zellkonverter)
library(DropletUtils)

# load arguments from the command line to this script
args = commandArgs(trailingOnly=TRUE)

# obtaining a H5AD file with scRNA-seq expressions.
example_h5ad <- system.file("extdata", args, package = "zellkonverter")
input_h5ad = readH5AD(example_h5ad)

# clearing expressions from duplicates
output_h5ad <- emptyDrops(assay(input_h5ad))

# saving of expressions cleared of duplicates
writeH5AD(output_h5ad, "output.h5ad")