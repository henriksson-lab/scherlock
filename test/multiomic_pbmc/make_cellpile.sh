java -jar ../../build/cellpile.jar build \
	reduced/10x_atac.cellpile \
	sizes.chromosome \
	reduced/atac_possorted_genome_bam.bam barcodes.tsv.gz


java -jar ../../build/cellpile.jar build \
	reduced/10x_gex.cellpile \
	sizes.chromosome \
	reduced/gex_possorted_genome_bam.bam barcodes.tsv.gz

