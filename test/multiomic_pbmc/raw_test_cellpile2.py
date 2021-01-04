### For cellpiles
import scherlock.CellPile as cpl
import scanpy as sc

# count tables
adata = sc.read_h5ad("reduced/10x_gex.h5ad")
adata.var_names_make_unique()
adata.obs["sampleid"] = "sample1"
adata.obs["leiden"] = "all"

# Initialize cellpile
cp = cpl.CellPile()

### Load the reference genome
cp.addGTF("reduced/ref/genes.gtf.gz")
cp.addPile("reduced/10x_atac.cellpile", pileName='sample1')  #Note that the pilename must correspond to sampleid


### Check our favourite genes, global pileup
v = cp.getView("CD55")

cp.pileup(v,
          cellBC=adata.obs["origbc"].tolist(),      # subset on cell BCs
          cellFile=adata.obs["sampleid"].tolist(),
          cellCluster=adata.obs["leiden"].tolist()  # Split into different tracks for clusters
         ).plot(save="out/cd55.svg")

