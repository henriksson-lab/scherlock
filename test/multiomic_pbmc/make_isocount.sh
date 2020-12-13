#Define gene features
java -jar ../../build/isocount.jar build reduced/ref/genes.gtf.gz reduced/featurefile.ff

#Build the count table (reduced)
java -jar ../../build/isocount.jar count reduced/isocount.gex reduced/featurefile.ff reduced/gex_possorted_genome_bam.bam barcodes.tsv.gz

#Build the count table (whole)
java -jar ../../build/isocount.jar count orig/isocount.gex reduced/featurefile.ff orig/gex/possorted_genome_bam.bam barcodes.tsv.gz
