
# cd to script directory so that relative paths work
# https://stackoverflow.com/questions/59895/how-can-i-get-the-source-directory-of-a-bash-script-from-within-the-script-itsel
sd="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$sd"


######################## reduced data ###########################################

java -jar ../../build/cellpile.jar build \
	'--single_cell' \
	reduced/10x_atac_single_cell.cellpile \
	sizes.chromosome \
	reduced/atac_possorted_genome_bam.bam \
	reduced/barcodes.tsv.gz

java -jar ../../build/cellpile.jar build \
	'--single_cell' \
	reduced/10x_gex_single_cell.cellpile \
	sizes.chromosome \
	reduced/gex_possorted_genome_bam.bam \
	reduced/barcodes.tsv.gz

java -jar ../../build/cellpile.jar build \
	'--bulk' \
	reduced/10x_atac_bulk.cellpile \
	sizes.chromosome \
	reduced/atac_possorted_genome_bam.bam \
	'bulk_sample'

java -jar ../../build/cellpile.jar build \
	'--bulk' \
	reduced/10x_gex_bulk.cellpile \
	sizes.chromosome \
	reduced/gex_possorted_genome_bam.bam \
	'bulk_sample'






# ######################## fullsize data ###########################################
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
