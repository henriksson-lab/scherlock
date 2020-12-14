all:

download_fullsize_test_data:
	echo
	echo 'NOTE:'
	echo 'This might take a while. You can also do'
	echo 'make download_reduced_test_data'
	echo 'to download a small test dataset'
	echo

	wget http://beagle.henlab.org/public/cellpile/test/orig.tar
	tar -zxvf orig.tar

download_reduced_test_data:
	wget http://beagle.henlab.org/public/cellpile/test/reduced.tar
	tar -zxvf reduced.tar


## Not functional/tested. Should we	even have this? Certainly not for users..
# upload_orig:
# 	mkdir -p orig
# 	rsync -r beagle.henlab.org:/data/public_http/public/cellpile/test/orig.tar .

# upload_reduced:
# 	mkdir -p reduced
# 	rsync -r beagle.henlab.org:/data/public_http/public/cellpile/test/reduced.tar .

# ToDo: Something nicer here..
chromsize:
	cut -f1,2 /big/anbjork/atac_pbmcs/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa.fai > sizes.chromosome

