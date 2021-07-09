### For cellpiles
import pathlib
import scherlock.cellpile as cpl
import scanpy as sc

# count tables
adata = sc.read_h5ad("reduced/10x_gex.h5ad")
adata.var_names_make_unique()
adata.obs["sampleid"] = "sample1"
adata.obs["leiden"] = "all"











##################### Single cell ######################### 

mode = 'single_cell'

for dataType in ['gex', 'atac']:

    outdir = "outs/" + mode + '/' + dataType

    # Initialize cellpile
    cp = cpl.cellpile()

    ### Load the reference genome
    cp.add_gtf("reduced/ref/genes.gtf.gz")
    # Note that the pilename must correspond to sample id
    cellpile_file = "reduced/10x_" + dataType + "_single_cell.cellpile"
    cp.add_cellpile(cellpile_file, cellpile_name='sample1')

    ### Check our favourite gene
    v = cp.get_view("CD55")

    # Plot and save
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

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
                 ).plot(save=(outdir + "/cd55__" + file_name_midfix + ".svg"))

    del cp  # Make sure the old cellpile object is not lingering













##################### Several tracks ######################### 

# Several tracks based on cell clusters is neccessarily single cell
mode = 'several_clusters'

aa = ['asdf' for elem in range(round(adata.shape[0]/2))]
bb = ['qwer' for elem in range(adata.shape[0] - round(adata.shape[0]/2))]
fake_clusters = aa + bb
assert(len(fake_clusters) == adata.shape[0])

for dataType in ['gex', 'atac']:

    outdir = "outs/" + mode + '/' + dataType

    # Initialize cellpile
    cp = cpl.cellpile()

    ### Load the reference genome
    cp.add_gtf("reduced/ref/genes.gtf.gz")
    # Note that the pilename must correspond to sample id
    cellpile_file = "reduced/10x_" + dataType + "_single_cell.cellpile"
    cp.add_cellpile(cellpile_file, cellpile_name='sample1')

    ### Check our favourite gene
    v = cp.get_view("CD55")

    # Plot and save
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

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
                  track_labels=fake_clusters,  # Split into different tracks for clusters
                  show_inbetweens=line[0],
                  individual_track_scaling=line[1],
                 ).plot(save=(outdir + "/cd55__" + file_name_midfix + ".svg"))

    del cp  # Make sure the old cellpile object is not lingering












##################### Bulk ######################### 


mode = 'bulk'

for dataType in ['gex', 'atac']:

    outdir = "outs/" + mode + '/' + dataType

    # Initialize cellpile
    cp = cpl.cellpile()

    ### Load the reference genome
    cp.add_gtf("reduced/ref/genes.gtf.gz")
    # Note that the pilename must correspond to sample id
    cellpile_file = "reduced/10x_" + dataType + "_bulk.cellpile"
    cp.add_cellpile(cellpile_file, cellpile_name='sample1')

    ### Check our favourite gene
    v = cp.get_view("CD55")

    # Plot and save
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

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
                  barcodes=["bulk_sample" for elem in range(adata.shape[0])],
                  cellpile_names=adata.obs["sampleid"].tolist(),
                  track_labels=adata.obs["leiden"].tolist(),
                  show_inbetweens=line[0],
                  individual_track_scaling=line[1],
                 ).plot(save=(outdir + "/cd55__" + file_name_midfix + ".svg"))

    del cp  # Make sure the old cellpile object is not lingering






