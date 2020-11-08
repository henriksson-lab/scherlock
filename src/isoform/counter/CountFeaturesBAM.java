package isoform.counter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * Count reads from a BAM file
 * 
 * @author Johan Henriksson
 *
 */
public class CountFeaturesBAM {
	//SAM format information:
	// https://genome.sph.umich.edu/wiki/SAM
	// https://www.biostars.org/p/96347/

	///////// Info for one read
	//-F INT   only include reads with none of the FLAGS in INT present [0]
	//read unmapped (0x4)
	//not primary alignment (0x100)
	//read fails platform/vendor quality checks (0x200)
	//read is PCR or optical duplicate (0x400)
	//supplementary alignment (0x800)
	//in total: 0xF04


	/**
	 * Dense count table. It will be compressed as we get through the sorted BAM-file and can rule out any addition to previous indices
	 */
	private int countTable[][];
	
	/**
	 * Compressed count table. Entries added line-by-line from the dense count table
	 */
	private int compressedCountTable[][];
	
	/**
	 * Map cell barcode -> array index
	 */
	private HashMap<String, Integer> mapBarcodeIndex=new HashMap<String, Integer>();

	private ArrayList<String> listBarcodes;
	private ArrayList<Feature> listFeatures;

	/**
	 * Compress feature at given index
	 */
	private void compressFeature(int i) {
		compressedCountTable[i]=compressFeature(countTable[i]);
		countTable[i]=null;
	}

	/**
	 * Compress an array into [index value, index value ...  ]
	 */
	private int[] compressFeature(int[] counts) {
		if(counts==null) {
			return null;
		} else {
			//Allocate an oversized array that can handle the worst case scenario
			int[] compressed=new int[counts.length*2];
			
			//Find non-zero counts
			int sparsei=0;
			for(int densei=0;densei<counts.length;densei++) {
				int c=counts[densei];
				if(c!=0) {
					compressed[sparsei]=densei;
					sparsei++;
					compressed[sparsei]=c;
					sparsei++;
				}
			}
			//Allocate a smaller array that exactly fits
			int[] smallerCompressed=new int[sparsei];
			System.arraycopy(compressed, 0, smallerCompressed, 0, sparsei);
			return smallerCompressed;
		}
	}
	
	/**
	 * Constructor
	 */
	public CountFeaturesBAM(ArrayList<Feature> listFeatures, ArrayList<String> listBarcodes) {
		this.listBarcodes=listBarcodes;
		this.listFeatures=listFeatures;
		
		//For quick-lookup of barcodes, if they should be considered
		for(int bi=0;bi<listBarcodes.size();bi++)
			mapBarcodeIndex.put(listBarcodes.get(bi), (Integer)bi);
		
		//Allocate the count table. Lazy dense format for simplicity
		System.out.println("Count matrix will have size: "+listFeatures.size()+ " * "+listBarcodes.size());
		countTable=new int[listFeatures.size()][];
		
		//Final place of storage. Features can be compressed as we have passed them completely in the count table
		compressedCountTable=new int[listFeatures.size()][];
	}
	
	/**
	 * Add one count. Allocate more dense space if needed
	 */
	private void count(int fi, int barcodeIndex) {
		if(countTable[fi]==null) {
			countTable[fi]=new int[listBarcodes.size()];
		}
		countTable[fi][barcodeIndex]++;
	}

	
	
	
	/**
	 * Get the next feature index that is at least not past the given block
	 */
	public int nextFeatureBeforeOrOverlap(String blockSource, int blockFrom, int searchFeatureListFrom) {
		//Move feature start search position to the next feature overlapping or at least not past the read
		for(;;) {
			//Stop if there are no more features to step through
			if(searchFeatureListFrom>=listFeatures.size()-1)
				break;

			//Check if we can step to the next feature or if we should remain
			Feature feature=listFeatures.get(searchFeatureListFrom+1);
			//System.out.println("Feature:\t"+feature.from+"\t"+feature.to);
			int comparisonSource=feature.source.compareTo(blockSource);  // could pre-split the chromosome to avoid these checks
			if(comparisonSource==0) {
				//Currently checking features on the same chromosome. Most common case.
				if(feature.from < blockFrom) {
					//Can still keep stepping
					searchFeatureListFrom++;
					continue;
				} else {
					//Cannot step further. We are done
					break;
				}
			} else if(comparisonSource<0) {
				//Not yet comparing the same chromosome. Keep scanning along the features until we get there. Uncommon case
				searchFeatureListFrom++;
			} else {
				//Now we are looking at the next chromosome, so past the block for certain.
				//This can happen if a read is mapped before any feature is defined on the chromosome.
				break;
			}
		}
		return searchFeatureListFrom;
	}
	
	
	
	
	
	/**
	 * Perform the counting from a BAM file. Can be called multiple times
	 */
	public void countReads(File fBAM) throws IOException {
		//Assume there is at least one feature to count
		if(listFeatures.isEmpty()) {
			System.out.println("No features to count");
			return;
		}

		//Open the BAM file and get to work
		final SamReader reader = SamReaderFactory.makeDefault().open(fBAM);
		
		
		//Statistics of how many records we kept
		int readRecords=0;
		int keptRecords=0;
		String currentSource=null;
		int currentPos=0;
		
		//Instead of keeping features in a treemap, which is O(log n) to look up,
		//we can exploit that the input file is sorted. By keeping a pointer to the
		//last feature used, we can do an O(1) lookup
		int searchFeatureListFrom=0;

		//This is to keep track of duplicates.
		//Approximate, as the same cDNA can be fragmented multiple times in the library prep
		String bcCellPreviousU=null;
		String bcCellPreviousC=null;
				
		//Loop through all SAM records
		for (final SAMRecord samRecord : reader) {
			readRecords++;
			if(readRecords%1000000 == 0){
				int prcDone=(int)(100.0*searchFeatureListFrom/(double)listFeatures.size());
				System.out.println("Progress: "+prcDone+"%\t  Kept/Read: "+keptRecords+"/"+readRecords+"\tCurrently@read: "+currentSource+"\t"+currentPos);
			}

			//Get UMI and BC for this read
			String bcCellCurrentUMI=(String)samRecord.getAttribute("UB");
			String bcCellCurrentCellBarcode=(String)samRecord.getAttribute("CB");
				
			//Check if this is a cell to count
			if(mapBarcodeIndex.keySet().contains(bcCellCurrentCellBarcode)) {
				//If the read has no UMI nor BC then ignore it
				if(bcCellCurrentUMI!=null && bcCellCurrentCellBarcode!=null) {
					//Check if duplicate read
					if(bcCellCurrentUMI.equals(bcCellPreviousU) && bcCellCurrentCellBarcode.equals(bcCellPreviousC)) {
						//Do nothing, just ignore - this is a duplicate read
					} else {
						//Remember for later
						bcCellPreviousU=bcCellCurrentUMI;
						bcCellPreviousC=bcCellCurrentCellBarcode;

						//Store position, for displaying progress
						currentSource=samRecord.getContig();
						currentPos=samRecord.getAlignmentStart();

						//Which cell is this?
						int barcodeIndex=mapBarcodeIndex.get(bcCellCurrentCellBarcode);
						
						//Move the search feature forward
						//System.out.println("Moving search initial position:");
						int nextSearchFeatureListFrom=nextFeatureBeforeOrOverlap(
								samRecord.getContig(), samRecord.getAlignmentStart(),
								searchFeatureListFrom);
						
						//Since the input is position sorted, we can compress all features up until now
						for(;searchFeatureListFrom<nextSearchFeatureListFrom;searchFeatureListFrom++) {
							compressFeature(searchFeatureListFrom);
						}
						

						/////////////////////////////////
						///////////////// Count the blocks
						/////////////////////////////////
						
						//A read may have been split into multiple blocks. 
						//Count these separately. Naive assumption that these are split over introns... is this correct?
						List<AlignmentBlock> listBlocks=samRecord.getAlignmentBlocks();
						for(int curAB=0;curAB<listBlocks.size();curAB++) {
							AlignmentBlock ab=listBlocks.get(curAB);
							
							String blockSource=samRecord.getContig();
							int blockFrom=ab.getReferenceStart();
							int blockTo=ab.getReferenceStart()+ab.getLength();

							//Find new suitable place to start searching from
							int fi=nextFeatureBeforeOrOverlap(samRecord.getContig(), ab.getReadStart(),searchFeatureListFrom);
							
							//Look up overlapping features. Continue searching from where we were last time.
							//This assumes that the FASTQ has been position sorted!
							//System.out.println("Fitting block");
							blocksearch: for(;fi<listFeatures.size();fi++) {
								Feature feature=listFeatures.get(fi);
								int comparisonSource=feature.source.compareTo(blockSource);  // could pre-split the chromosome to avoid these checks
								if(comparisonSource==0) {
									//Currently checking features on the same chromosome. Most common case.
									
									/*System.out.println("AB!\t"+blockSource+"\t"+blockFrom+"\t"+blockTo+"\t"+curAB+"\tumi "+bcCellCurrentUMI+"\tbc "+bcCellCurrentCellBarcode);
									System.out.println("Feature:\t"+feature.from+"\t"+feature.to);
									System.out.println();*/
									
									//Now continue to check if there is an overlap position-wise
									if(feature.to<blockFrom) {
										//Feature is before block. We have not reached the region yet. Continue the search
									} else if(feature.from>blockTo) {
										//Now we are past the feature, but on the same chromosome. Can stop looking
										break blocksearch;
									} else {
										//This feature overlaps, so count it. But continue the loop as multiple features might overlap the block
										count(fi,barcodeIndex);
										keptRecords++;
									}
								} else if(comparisonSource<0) {
									//Not yet comparing the same chromosome. Keep scanning along the features until we get there. Uncommon case
								} else {
									//Now we are looking at the next chromosome, so past the block for certain. Can stop looking
									break blocksearch;
								}
							}
						}
					}
				}
			}
		}

		//Need to ensure the final/current block is compressed and chucked away as well
		if(countTable[searchFeatureListFrom]!=null) {
			compressFeature(searchFeatureListFrom);
		}
		
		System.out.println("Kept/Read: "+keptRecords+"/"+readRecords+" --- "+(int)(100*keptRecords/(double)readRecords)+"%");

		reader.close();
	}
	

	/**
	 * Write the output in 10x format
	 */
	public void writeMatrix(File fCountDir) throws IOException {
		
		//Write the barcodes
		PrintWriter pwBarcodes=new PrintWriter(new GZIPOutputStream(new FileOutputStream(new File(fCountDir,"barcodes.tsv.gz"))));
		for(String bc:listBarcodes) {
			pwBarcodes.println(bc);
		}
		pwBarcodes.close();

		//Write the feature names
		PrintWriter pwFeatures=new PrintWriter(new GZIPOutputStream(new FileOutputStream(new File(fCountDir,"features.tsv.gz"))));
		for(Feature f:listFeatures) {
			pwFeatures.println(f.featureName+"\t"+f.featureName+"\tisoform");
		}
		pwFeatures.close();
		
		//Count non-zero counts; needed for matrix market
		long numNonZero=0;
		for(int fi=0;fi<listFeatures.size();fi++) {
			int[] countsForF=compressedCountTable[fi];
			if(countsForF!=null) {
				numNonZero+=countsForF.length/2;
			}
		}

		
		//Write the count matrix
		PrintWriter pwMatrix=new PrintWriter(new GZIPOutputStream(new FileOutputStream(new File(fCountDir,"matrix.mtx.gz"))));
		pwMatrix.println("%%MatrixMarket matrix coordinate integer general");
		pwMatrix.println("%metadata_json: {\"format_version\": 2, \"software_version\": \"3.1.0\"}");  //well...
		pwMatrix.println(""+listFeatures.size()+" "+listBarcodes.size()+" "+numNonZero);
		for(int fi=0;fi<listFeatures.size();fi++) {
			//Write a compressed feature count
			int[] countsForFeature=compressedCountTable[fi];
			if(countsForFeature!=null) {
				for(int i=0;i<countsForFeature.length;i+=2) {
					int bcIndex=countsForFeature[i];
					int count=countsForFeature[i+1];
					pwMatrix.println(""+(fi+1)+" "+(bcIndex+1)+" "+count);
					//Note that matrix market indices are 1-based
				}
			}
		}
		pwMatrix.close();
	}
	
	
	
	
	
}
