import scanpy as sc
import pandas as pd

list_bc = pd.read_csv("barcodes.tsv.gz",header=None)[0].tolist()

###############################################################
### Subset 10x GEX matrix to reduced set of barcodes ##########
###############################################################

one_sample = sc.read("orig/gex/matrix.mtx.h5ad")
sub_sample = one_sample[[x in list_bc for x in one_sample.obs_names],:]
sub_sample.write_h5ad("reduced/10x_gex.h5ad")

###############################################################
### Subset 10x GEX matrix to reduced set of barcodes ##########
###############################################################

one_sample = sc.read_10x_h5("orig/atac/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
sub_sample = one_sample[[x in list_bc for x in one_sample.obs_names],:]
sub_sample.write_h5ad("reduced/10x_atac.h5ad")
