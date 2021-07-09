### For cellpiles
import pathlib
import scherlock.cellpile as cpl
import scanpy as sc

# count tables
adata = sc.read_h5ad("reduced/10x_gex.h5ad")
adata.var_names_make_unique()
adata.obs["sampleid"] = "sample1"
adata.obs["leiden"] = "all"

# Initialize cellpile
cp = cpl.cellpile()

### Load the reference genome
cp.add_gtf("reduced/ref/genes.gtf.gz")
# Note that the pilename must correspond to sample id
cp.add_cellpile("reduced/10x_atac.cellpile", cellpile_name='sample1')  

### Check our favourite gene
v = cp.get_view("CD55")

# Plot and save
pathlib.Path("outs/").mkdir(parents=True, exist_ok=True)

testing_order = [[False, False],
                 [False, True ],
                 [True,  False],
                 [True,  True ]]


for line in testing_order:
    testing_args = {"show_inbetweens" : line[0],
                    "individual_track_scaling" : line[1]}
    file_name_midfix = "".join([key + '_' + str(value).ljust(5, '_') + '___' 
        for key, value in testing_args.items()]).strip('_')
    cp.pileup(v,
              barcodes=adata.obs["origbc"].tolist(),      # subset on cell BCs
              cellpile_names=adata.obs["sampleid"].tolist(),
              track_labels=adata.obs["leiden"].tolist(),  # Split into different tracks for clusters
              show_inbetweens=line[0],
              individual_track_scaling=line[1],
             ).plot(save="outs/cd55__" + file_name_midfix + ".svg")

