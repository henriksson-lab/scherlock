import scanpy as sc
import pandas as pd
import scrublet as scr
import scherlock

######################### Read the data

adata = sc.read_h5ad("orig/gex/matrix.mtx.h5ad")

# but reduced ... 500



######################### Cell tables


### Calculate absolute cell counts
print(cellcounts(adata, ["leiden"]))

### Calculate percentage cells of the whole
print(cellcounts(adata, ["leiden"], percentage=True))

### Calculate percentage cells, within each group
#cellcounts(adata, ["genotype", "treatment"], percentage=True, groupby=["genotype"])



