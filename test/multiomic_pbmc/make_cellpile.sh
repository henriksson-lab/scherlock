######################## reduced ###########################################
java -jar ../../build/cellpile.jar build \
	'--single_cell' \
	reduced/10x_atac.cellpile \
	sizes.chromosome \
	reduced/atac_possorted_genome_bam.bam \
	reduced/barcodes.tsv.gz


java -jar ../../build/cellpile.jar build \
	'--single_cell' \
	reduced/10x_gex.cellpile \
	sizes.chromosome \
	reduced/gex_possorted_genome_bam.bam \
	reduced/barcodes.tsv.gz

# ######################## original ###########################################
# java -jar ../../build/cellpile.jar build \
# 	fullsize/10x_atac.cellpile \
# 	sizes.chromosome \
# 	fullsize/atac/possorted_genome_bam.bam fullsize/barcodes.tsv.gz


# java -jar ../../build/cellpile.jar build \
# 	fullsize/10x_gex.cellpile \
# 	sizes.chromosome \
# 	fullsize/gex/possorted_genome_bam.bam fullsize/barcodes.tsv.gz
# # # Old one. Shouldnt be reduced, right?
# # java -jar ../../build/cellpile.jar build \
# # 	reduced/10x_gex.cellpile \
# # 	sizes.chromosome \
# # 	fullsize/gex/possorted_genome_bam.bam fullsize/barcodes.tsv.gz
