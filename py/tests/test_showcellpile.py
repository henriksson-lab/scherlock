import scherlock


cp = scherlock.cp.CellPile()
cp.addPile("/home/mahogny/temp/cellpile")
cp.loadGTF("/home/mahogny/umeå/project/isoform/refgenome/Homo_sapiens.GRCh38.101.chr.gff3.gz")


#### installation:
#pip install py4j




grange = cp.getView("CD55")
print(grange)


#cd55
out = cp.pileup(grange)
out.plot(save="/home/mahogny/temp.svg")


