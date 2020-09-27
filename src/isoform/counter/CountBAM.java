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
public class CountBAM {
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


	private int countTable[][];
	
	private int compressedCountTable[][];
	
	
	private HashMap<String, Integer> mapBarcodeIndex=new HashMap<String, Integer>();

	private ArrayList<String> listBarcodes;
	private ArrayList<Feature> listFeatures;
	
	private int dups=0;
	private int kept=0;

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
					compressed[sparsei]=densei+1; //MTX market starts from 1
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
	public CountBAM(ArrayList<Feature> listFeatures, ArrayList<String> listBarcodes) {
		this.listBarcodes=listBarcodes;
		this.listFeatures=listFeatures;
		
		//For quick-lookup of barcodes, if they should be considered
		//HashSet<String> setBarcodes=new HashSet<String>(listBarcodes);
		for(int bi=0;bi<listBarcodes.size();bi++)
			mapBarcodeIndex.put(listBarcodes.get(bi), (Integer)bi);
		
		
		//Allocate the count table. Lazy dense format for simplicity
		System.out.println("Count matrix will have size: "+listFeatures.size()+ " * "+listBarcodes.size());
		countTable=new int[listFeatures.size()][];//[listBarcodes.size()];
		
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
		kept++;
	}

	
	/**
	 * Perform the counting from a BAM file. Can be called multiple times
	 */
	public void countReads(File fBAM) throws IOException {
		final SamReader reader = SamReaderFactory.makeDefault().open(fBAM);
		
		int readRec=0;
		int searchFeatureFrom=0;

		SAMRecord previousRecord=null;
		String bcCellPreviousU=null;
		String bcCellPreviousC=null;
				
		for (final SAMRecord samRecord : reader) {
			readRec++;
			if(readRec%1000000 == 0){
				int prcDone=(int)(100.0*searchFeatureFrom/(double)listFeatures.size());
				System.out.println("Progress: "+prcDone+"%\t  Kept/Read"+kept+"/"+readRec);
			}

			if(previousRecord!=null) {
				String bcCellCurrentUMI=(String)samRecord.getAttribute("UB");
				String bcCellCurrentCellBarcode=(String)samRecord.getAttribute("CB");
				
				//Check if this is a cell to count
				if(mapBarcodeIndex.keySet().contains(bcCellCurrentCellBarcode)) {
					//If the read has no UMI nor BC then ignore it
					if(bcCellCurrentUMI!=null && bcCellCurrentCellBarcode!=null) {
						//Check if duplicate read -- simplistic method but maybe good enough?
						if(bcCellCurrentUMI.equals(bcCellPreviousU) && bcCellCurrentCellBarcode.equals(bcCellPreviousC)) {
							dups++;
						} else {
							//Remember for later
							bcCellPreviousU=bcCellCurrentUMI;
							bcCellPreviousC=bcCellCurrentCellBarcode;

							//Which cell is this?
							int barcodeIndex=mapBarcodeIndex.get(bcCellCurrentCellBarcode);
							
							//A read may have been split into multiple blocks. 
							//Count these separately. Naive assumption that these are split over introns... is this correct?
							List<AlignmentBlock> listBlocks=samRecord.getAlignmentBlocks();
							for(int curAB=0;curAB<listBlocks.size();curAB++) {
								AlignmentBlock ab=listBlocks.get(curAB);
								
								String blockSource=samRecord.getContig();
								int blockFrom=ab.getReferenceStart();
								int blockTo=ab.getReferenceStart()+ab.getLength();
								
								//Look up overlapping feature
								for(int fi=searchFeatureFrom;fi<listFeatures.size();fi++) {
									Feature feature=listFeatures.get(fi);
									int compSource=feature.source.compareTo(blockSource);  // could pre-split the chromosome to avoid these checks
									if(compSource==0) {
										//Currently checking features on the same chromosome. Most common case.
										
										//Now continue to check if there is an overlap position-wise
										if(feature.to<blockFrom) {
											//Feature is before block. We have not reached the region yet.
											
											//DANGER: are alignment blocks in left-to-right order? VERIFY!!!!!!!!!!!!!!!!! OTHERWISE SORT THEM FIRST!
											
											//Since the input is position sorted, if we reach this code then this feature will never
											//be counted again. Thus, we can compress all features up until now
											for(;searchFeatureFrom<=fi;searchFeatureFrom++) {
												compressedCountTable[searchFeatureFrom]=compressFeature(countTable[searchFeatureFrom]);
												countTable[searchFeatureFrom]=null;
											}
											
											continue;
										} else if(feature.from>blockTo) {
											//Now we are past the feature. Can stop looking
											break;
										} else {
											//This feature overlaps, so count it
											count(fi,barcodeIndex);
											/*
											System.out.println("---overlap!");
											System.out.println(blockSource+":"+blockFrom+"-"+blockTo);
											System.out.println(feature.source+":"+feature.from+"-"+feature.to);
											System.out.println(feature.gene);
											System.out.println("---");
											System.out.println(blockSource+"\t"+blockFrom+"\t"+blockTo+"\t"+"theread");
											System.out.println(feature.source+"\t"+feature.from+"\t"+feature.to+"\t"+feature.featureName);
											System.exit(0);*/
											
											//Here we do *not* Continue;, as there might be more overlapping features 
										}
									} else if(compSource<0) {
										//Not yet comparing the same chromosome. Keep scanning along the features until we get there. Uncommon case
									} else {
										//Now we are looking at the next chromosome, so past the block for certain. Can stop looking
										break;
									}
								}
							}
						}
					}
				}
			} else {
				previousRecord=samRecord;
			}
		}
		System.out.println("Kept: "+kept);
		System.out.println("Duplicates: "+dups);

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
		
		//Write the count matrix
		PrintWriter pwMatrix=new PrintWriter(new GZIPOutputStream(new FileOutputStream(new File(fCountDir,"matrix.mtx.gz"))));
		pwMatrix.println("%%MatrixMarket matrix coordinate integer general");
		pwMatrix.println("%metadata_json: {\"format_version\": 2, \"software_version\": \"3.1.0\"}");  //well...
		pwMatrix.println(""+listFeatures.size()+" "+listBarcodes.size()+" "+listBarcodes.size());
		for(int fi=0;fi<listFeatures.size();fi++) {
			
			
			//Write a compressed feature count
			int[] countsForF=compressedCountTable[fi];
			if(countsForF!=null) {
				for(int i=0;i<countsForF.length;i+=2) {
					int bi=countsForF[i];
					int c=countsForF[i+1];
					pwMatrix.println(""+(fi+1)+" "+(bi+1)+" "+c);
				}
			}
			
			//Below is for the dense matrix. But now it is compressed while computing
			/*
			int[] countsForF=countTable[fi];
			for(int bi=0;bi<listBarcodes.size();bi++) {
				int c=countsForF[bi];
				if(c!=0) {
					pwMatrix.println(""+(fi+1)+" "+(bi+1)+" "+c);				
				}
			}
			*/
		}
		pwMatrix.close();
	}
	
	
	
	
	
}
