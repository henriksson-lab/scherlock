import scanpy as sc
import pandas as pd
import scrublet as scr
import scherlock
from random import random

######################### Read the data

adata = sc.read_h5ad("reduced/10x_gex_500.h5ad")

#invent column values
adata.obs["treatment"] = [x<250 for x in range(0,adata.shape[0])]
adata.obs["genotype"] = ["ko" if random()<0.5 else "wt" for x in range(0,adata.shape[0])]


######################### Cell tables


### Calculate absolute cell counts
print(scherlock.stats.cellcounts(adata, ["celltype"]))

### Calculate percentage cells of the whole
print(scherlock.stats.cellcounts(adata, ["celltype"], percentage=True))

### Calculate percentage cells, within each group
print(scherlock.stats.cellcounts(adata, ["genotype", "treatment"], percentage=True, groupby=["genotype"]))




