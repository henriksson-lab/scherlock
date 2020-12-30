

set -e
trap "kill 0" SIGINT


echo 'Started:'
date
echo

java -jar ../../build/cellpile.jar randomsubset \
	reduced/barcodes.tsv.gz \
	fullsize/gex/possorted_genome_bam.bam \
	reduced/gex_possorted_genome_bam.bam  0.5 0.0005 \
    &> reduce_bam.log_gex &

java -jar ../../build/cellpile.jar randomsubset \
	reduced/barcodes.tsv.gz \
	fullsize/atac/possorted_genome_bam.bam \
	reduced/atac_possorted_genome_bam.bam  0.5 0.0005 \
    &> reduce_bam.log_atac &

wait


echo 'Done :)'
date
echo


