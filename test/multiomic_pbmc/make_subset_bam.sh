java -jar ../../build/cellpile.jar randomsubset \
	gex/barcodes.tsv.gz \
	/big/anbjork/atac_pbmcs/data/gex/fake_cellranger_sample_output_directory/possorted_genome_bam.bam \
	reduced/gex_possorted_genome_bam.bam  0.5 0.0005

java -jar ../../build/cellpile.jar randomsubset \
	gex/barcodes.tsv.gz \
	/big/anbjork/atac_pbmcs/data/atac/download/possorted_genome_bam.bam \
	reduced/atac_possorted_genome_bam.bam  0.5 0.0005

########################


