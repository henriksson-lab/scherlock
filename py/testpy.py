import python_interface


cp = CellPile("/home/mahogny/temp/cellpile")

#should support multiple cellpiles somehow

print(cp.getSequences())

#### installation:
#pip install py4j


#//cd55: "1",207321376,207360966

#Some command to get a range around a gene

cp.loadGTF("/home/mahogny/umeå/project/isoform/refgenome/Homo_sapiens.GRCh38.101.chr.gff3.gz")

grange = cp.getView("CD55")
#grange = ("1", 207321376,207360966)
print(grange)


#cellgroups = ...   #quick function to make one single group

#need features whenever multiple functons have been used

#cd55
cp.pileup(grange, cellgroups)