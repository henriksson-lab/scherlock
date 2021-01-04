import numpy as np
import scipy
import pandas as pd
import anndata
import scanpy as sc
import scipy.io
#import matplotlib
#from matplotlib.backends.backend_pdf import PdfPages
#import matplotlib.pyplot as plt
import os


sc.settings.verbosity = 3


#########################
######################### Read the data (pre-annotated) ... or will be
#########################

adata = sc.read_h5ad("reduced/10x_gex.h5ad")

adata.var_names_make_unique()
adata.obs["sampleid"] = "sample1"

# annotate the group of mitochondrial genes as 'mt'
adata.var['mt'] = adata.var["gene_names"].str.startswith('MT-')  
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)

# Normalize the counts
sc.pp.normalize_total(adata, target_sum=1e6)    

# Move to log counts
sc.pp.log1p(adata) 

### Store what we have now, somewhat raw (CPM-normalized and pseudo-log)
adata.raw=adata




def findHighlyVariable(adata):
    adata2 = sc.AnnData(X=adata.raw.X, 
                        var=adata.raw.var,
                        obs = adata.obs)
    mitochondrial_genes = adata2.var_names[adata2.var["gene_names"].str.startswith('MT-')].tolist()

    #Figure out the variable genes for the full dataset
    sc.pp.highly_variable_genes(adata2, min_mean=0.5, max_mean=20, min_disp=0.3, batch_key = 'sampleid')
    sc.pl.highly_variable_genes(adata2)

    # Exclude the mitochondrial genes
    for i,g in enumerate(adata2.var.highly_variable.index.values.tolist()):
        if g in mitochondrial_genes:
            adata2.var.highly_variable[g] = False

    sc.pl.highly_variable_genes(adata2)

    ## Set which genes to consider highly variable
    var_genes = adata2.var.highly_variable[adata2.var.highly_variable == True].copy().index
    adata.var.highly_variable=adata2.var.highly_variable

    return(var_genes)

var_genes=findHighlyVariable(adata)




#Calculate neighbour graph
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.umap(adata, n_components=3)


### Perform the clustering
sc.tl.leiden(adata, resolution=0.3)


### For cellpiles
import scherlock.CellPile as cpl

# Initialize cellpile
cp = cpl.CellPile()

### Load the reference genome
cp.addGTF("reduced/ref/genes.gtf.gz")
cp.addPile("reduced/10x_atac.cellpile", pileName='sample1')  #Note that the pilename must correspond to sampleid


### Check our favourite genes, global pileup
v = cp.getView("CD55")

cp.pileup(v,
          cellBC=adata.obs["origbc"].tolist(),  # subset on cell BCs
          cellFile=adata.obs["sampleid"].tolist(), 
          cellCluster=adata.obs["leiden"].tolist()  # Split into different tracks for clusters
         ).plot(save="out/basic_cellpile.svg")

