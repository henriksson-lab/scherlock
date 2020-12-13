######################## reduced ###########################################
java -jar ../../build/cellpile.jar build \
	reduced/10x_atac.cellpile \
	sizes.chromosome \
	reduced/atac_possorted_genome_bam.bam reduced/barcodes.tsv.gz


java -jar ../../build/cellpile.jar build \
	reduced/10x_gex.cellpile \
	sizes.chromosome \
	reduced/gex_possorted_genome_bam.bam reduced/barcodes.tsv.gz

######################## original ###########################################
java -jar ../../build/cellpile.jar build \
	orig/10x_atac.cellpile \
	sizes.chromosome \
	orig/atac/possorted_genome_bam.bam orig/barcodes.tsv.gz


java -jar ../../build/cellpile.jar build \
	reduced/10x_gex.cellpile \
	sizes.chromosome \
	orig/gex/possorted_genome_bam.bam orig/barcodes.tsv.gz
