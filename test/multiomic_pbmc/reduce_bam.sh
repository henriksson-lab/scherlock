java -jar ../../build/cellpile.jar randomsubset \
	reduced/barcodes.tsv.gz \
	orig/gex/possorted_genome_bam.bam \
	reduced/gex_possorted_genome_bam.bam  0.5 0.0005

java -jar ../../build/cellpile.jar randomsubset \
	reduced/barcodes.tsv.gz \
	orig/atac/possorted_genome_bam.bam \
	reduced/atac_possorted_genome_bam.bam  0.5 0.0005

########################


