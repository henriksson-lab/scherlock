package isoform.cellpile;
import java.lang.Integer;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Set;
import java.util.HashMap;
import java.util.List;
import java.util.StringTokenizer;
import java.util.TreeMap;
import cern.colt.list.IntArrayList;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import isoform.util.LogUtil;
import isoform.util.PileUtil;

/**
 * Interface for CellPile files. Able to generate from BAM-files
 * 
 * Consider swapping to a buffered implementation for the read interface (need to verify if this improves speed)
 * 
 * @author Johan Henriksson and Anton Björk
 *
 */
public class CellPileFile {
	
	private RandomAccessFile raf;
	
	//All of this is in the header
	private int chunkSize=10000; //10kb
	private TreeMap<String, long[]> mapChunkStarts=new TreeMap<String, long[]>();
	public HashMap<String, Integer> mapBarcodeIndex=new HashMap<String, Integer>();
	private ArrayList<String> listBarcodes;
	
	/**
	 * Hidden constructor, cannot be called from the outside.
	 * Use the read/write functions instead
	 */
	private CellPileFile() {}
	

	
	/////////////////////////////////////////////////////////////////////////////////////
	///////////////////// Write CellPiles   /////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////

	/////////////////////
	//For building: Current set of regions that we are collecting
	private TreeMap<Integer, IntArrayList> mapCellAlignmentBlocks=new TreeMap<Integer, IntArrayList>();
	private TreeMap<Integer, IntArrayList> mapCellInbetweens=new TreeMap<Integer, IntArrayList>();
	private TreeMap<Integer, IntArrayList> mapCellLeftoverAlignmentBlocks=new TreeMap<Integer, IntArrayList>();
	private TreeMap<Integer, IntArrayList> mapCellLeftoverInbetweens=new TreeMap<Integer, IntArrayList>();

	private String currentSeq;
	private int currentChunk;

	
	/**
	 * Write a string to file
	 */
	public void writeString(String s) throws IOException {
		raf.writeShort((short)s.length());
		for(char c:s.toCharArray())
			raf.writeChar(c);
	}
	
	/**
	 * Read chromosome sizes and calculate chunk sizes
	 */
	// mapChunkStarts contains the starting points of chunks within the
	// different chromosomes. Note that ech chromosome has several chunks
	private void readChromosomeSizes(File f) throws IOException {
		BufferedReader br=new BufferedReader(new FileReader(f));
		String line;
		while((line=br.readLine())!=null) {
			StringTokenizer stok=new StringTokenizer(line,"\t");
			String seq=stok.nextToken();
			int len=Integer.parseInt(stok.nextToken());
			
			int numChunk=(len/chunkSize)+1;
			mapChunkStarts.put(seq,new long[numChunk]);
		}
		br.close();
	}
	
	/**
	 * Take a BAM-file, write a cellpile file corresponding to it
	 */
	public static CellPileFile writeFile(
			File fCellpile, 
			File fChromSizes,
			File fBAM,
			String bamType,
			ArrayList<String> listBarcodes) throws IOException {
		
		//Delete file or it will keep increasing in size
		if(fCellpile.exists())
			fCellpile.delete();
		
		CellPileFile p=new CellPileFile();
		p.readChromosomeSizes(fChromSizes);
		p.raf=new RandomAccessFile(fCellpile, "rw");
		p.listBarcodes=listBarcodes;
		
		//Write header - this gets the position right for later writes
		p.writeHeader();

		p.countReads(fBAM, bamType);
		
		//Update header with pointers
		p.writeHeader();

		return p;
	}
	
	/**
	 * For quick-lookup of barcodes, if they should be considered
	 */
	private void prepareBarcodeIndexLookup() {
		for(int bi=0;bi<listBarcodes.size();bi++)
			mapBarcodeIndex.put(listBarcodes.get(bi), (Integer)bi);
	}
	
	/**
	 * Write the header to the file
	 */
	private void writeHeader() throws IOException {
		//Go to the beginning of the file
		raf.seek(0);
		
		//Write the magic string and version
		byte version=2;
		raf.write(new byte[] {'C','E','L','L','P','I','L','E',version});
		
		raf.writeInt(chunkSize);
		
		//Write the names of cells
		prepareBarcodeIndexLookup();
		int numBarcodes=mapBarcodeIndex.size();
		raf.writeInt((int)numBarcodes);
		for(String s:listBarcodes) {
			writeString(s);
		}
		
		//Write all the chunk start positions
		int numSequences=mapChunkStarts.size();
		raf.writeInt(numSequences);
		for(String seqname:mapChunkStarts.keySet()) {
			long[] chunkstart=mapChunkStarts.get(seqname);
			writeString(seqname);
			raf.writeInt(chunkstart.length);
			for(long v:chunkstart)
				raf.writeLong(v);
		}
	}
	
	

	/**
	 * Add a region for a given barcode/cell
	 */
	private void addRegion(TreeMap<Integer, IntArrayList> regionType, 
							int bcid, int from, int to) {

		IntArrayList regs=regionType.get(bcid);
		if(regs==null) {
			regionType.put(bcid, regs=new IntArrayList(1000));
		}
		regs.add(from);
		regs.add(to);
	}
	

	/**
	 * Save current chunk
	 */
	private void saveCurrentChunk() throws IOException {

		int currentChunkFrom = currentChunk * chunkSize;
		int currentChunkTo = currentChunkFrom + chunkSize;

		addFromLeftovers(mapCellAlignmentBlocks, mapCellLeftoverAlignmentBlocks,
						 currentChunkFrom, currentChunkTo);
		addFromLeftovers(mapCellInbetweens, mapCellLeftoverInbetweens,
						 currentChunkFrom, currentChunkTo);

        // Store chunk to file
		// Only store chunks with data. Pointer is default 0=missing
		if(!mapCellAlignmentBlocks.isEmpty()) {
			//Ensure that this chunk has been allocated in the table; otherwise ignore it
			long[] ptrs=mapChunkStarts.get(currentSeq);
			if(ptrs!=null) {
				if(currentChunk<ptrs.length) {
					//Store pointer to where this chunk starts in the file
					ptrs[currentChunk]=raf.getFilePointer();

					// Store the chunk. This includes both the alignment blocks
					// and the inbetweens. No need for separate pointers since
					// they will typically be read together anyways.
					int numCellsAlignmentBlocks=mapCellAlignmentBlocks.size();
					int numCellsInbetweens=mapCellInbetweens.size();
					
					raf.writeInt(numCellsAlignmentBlocks);
					saveRegionTypeToFile(mapCellAlignmentBlocks);
					
					raf.writeInt(numCellsInbetweens);
					saveRegionTypeToFile(mapCellInbetweens);
					
					//System.out.println("Stored chunk, "+currentChunk+ " #cells "+mapCellAlignmentBlocks.size());
				} else {
					System.out.println("Warning, chunk out of bounds and ignored: "+currentSeq+":"+currentChunk);
				}
			} else {
				System.out.println("Warning, missing sequence ignored: "+currentSeq);
			}
		} else {
			System.out.println("No cells to write for chunk "+currentChunk);
		}
	}


	// Adds to current chunk from leftovers from previous chunks.
	// Just to avoid code duplication for the two region types in saveCurrentChunk()
	private void addFromLeftovers(TreeMap<Integer, IntArrayList> regionType,
								  TreeMap<Integer, IntArrayList> leftoverType,
								  int currentChunkFrom,
								  int currentChunkTo) throws IOException {
		Set<Integer> cellBarcodes = leftoverType.keySet();
        for(int bcid: cellBarcodes){
        	IntArrayList aa=leftoverType.get(bcid);
			for(int ii=0; ii<aa.size(); ii=ii+2) {
				int from = aa.get(ii);
				int to = aa.get(ii+1);
				if ((currentChunkFrom <= from) && (from < currentChunkTo)) {
					addRegion(regionType, bcid, from, to);
					// Remove entries to avoid leftovers growing large
					aa.remove(ii+1);
					aa.remove(ii);
				}
			}
        }
	}

	// Saves one region type (stored in one TreeMap) to file.
	// Just to avoid code duplication for the two region types in saveCurrentChunk()
	private void saveRegionTypeToFile(TreeMap<Integer, IntArrayList> regionType) throws IOException {
		for(int cellBarcode:regionType.keySet()) {
			raf.writeInt(cellBarcode);
			IntArrayList ia=regionType.get(cellBarcode);
			//Note, this is not a copy of the array, and it is the wrong size
			int[] list=ia.elements(); 
			int numRegion2=ia.size();
			raf.writeShort(numRegion2/2);
			for(int i=0;i<numRegion2;i++) {
				raf.writeInt(list[i]);
			}
		}
	}





	/**
	 * Save current chunk and move to the next position
	 */
	private void saveAndSetCurrentChunk(String newSeq, int newChunkNum) throws IOException {
		
		saveCurrentChunk();
		
		//Reset chunk store so we can fill with the next region
		mapCellAlignmentBlocks.clear();
		mapCellInbetweens.clear();
		
		// Clear the lefover buckets when switching chromosome since
		// mapCellAlignmentBlocks and mapCellInbetweens are not 
		// chromosome aware, and leftover buckets are used to append them
		// This introduces a minor loss of correctness when alignmentblocks
		// from one read span several chromosomes.
		// This should be very unusual though.
		if (currentSeq != newSeq) {
			mapCellLeftoverAlignmentBlocks.clear();
			mapCellLeftoverInbetweens.clear();
		}

		currentSeq=newSeq;
		currentChunk=newChunkNum;
	}

	/**
	 * Perform the counting from a BAM file. 
	 Can be called multiple times.  
	 Why call several times though? //AB
	 */
	public void countReads(File fBAM, String bamType) throws IOException {
		//Open BAM file
		final SamReader reader = SamReaderFactory.makeDefault().open(fBAM);

		//How many chunks are there?
		int totalChunks=0;
		for(String seq:mapChunkStarts.keySet())
			totalChunks+=mapChunkStarts.get(seq).length;

		//Set initial chunk position to nowhere
		currentSeq="";
		currentChunk=0;

		//To know progress, status variables
		int readRecords=0;
		int keptRecords=0;
		int skippedWrongBC=0;
		int skippedBadUMI=0;
		int skippedDup=0;
		int skippedNoAlignmentBlocks=0;

		//This is to keep track of duplicates.
		//Approximate, as the same cDNA can be fragmented multiple times in the library prep
		//    Assumes that positions are identical of the duplicates so that they will
		//    be directly after each other in the *position sorted* input BAM file
		String bcCellPreviousUMI="";
		String bcCellPreviousCB="";


		if (bamType.equals("--single_cell")) {
	
			//Loop through all SAM records
			for (final SAMRecord samRecord : reader) {
				
				//Update user about progress
				readRecords++;
				if(readRecords%1000000 == 0){
					//Calculate progress
					int passedChunks=0;
					for(String seq:mapChunkStarts.keySet()) {
						if(seq.equals(currentSeq))
							break;
						passedChunks+=mapChunkStarts.get(seq).length;
					}
					passedChunks+=currentChunk;
					int prc=100*passedChunks/totalChunks;
					
					LogUtil.formatColumns(System.out, 25,
							prc+"%",
							"Kept/Read: "+keptRecords+"/"+readRecords,
							"@sequence: "+currentSeq,
							"WrongBC: "+skippedWrongBC,
							"badUMI: "+skippedBadUMI,
							"SkipDup: "+skippedDup,
							"noAlignBlocks: "+skippedNoAlignmentBlocks);
				}
					
				//Get UMI and BC for this read
				String bcCellCurrentUMI=(String)samRecord.getAttribute("UB");
				String bcCellCurrentCellBarcode=(String)samRecord.getAttribute("CB");
				
				//Check if this is a cell to count
				if(mapBarcodeIndex.keySet().contains(bcCellCurrentCellBarcode)) {
					//If the read has no UMI nor BC then ignore it
					if(bcCellCurrentCellBarcode!=null) {
						//Check if duplicate read, if UMI present; ATAC, dedup by coordinate?
						if(bcCellCurrentUMI!=null && 
							bcCellCurrentUMI.equals(bcCellPreviousUMI) && 
							bcCellCurrentCellBarcode.equals(bcCellPreviousCB)) {
							skippedDup++;
						} else {
							//Remember for later
							bcCellPreviousUMI=bcCellCurrentUMI;
							bcCellPreviousCB=bcCellCurrentCellBarcode;

							//Which cell is this?
							int barcodeIndex=mapBarcodeIndex.get(bcCellCurrentCellBarcode);

							//A read may have been split into multiple blocks.
							List<AlignmentBlock> listBlocks=samRecord.getAlignmentBlocks();
							if (listBlocks.size() == 0) { 
								skippedNoAlignmentBlocks++;
								continue;  	// Wouldn't it be a good idea to use
										   	// more continue statements instead
											// of all the current nested if else?
											// if .. continue combo doesn't
											// introduce deeper nesting.
											// //AB
							}

							String blockSource=samRecord.getContig();

							// If read starts outside chunk, change chunk,
							// since BAM is position sorted on start position
							// of reads.
							// "leftmost coordinates"
							// http://www.htslib.org/doc/samtools-sort.html
							int aa=listBlocks.get(0).getReferenceStart();
							int shouldBeInChunk=aa/chunkSize;
							if((!currentSeq.equals(blockSource) || 
							   currentChunk!=shouldBeInChunk)) {
								saveAndSetCurrentChunk(blockSource, shouldBeInChunk);
							}

							// Add the alignment blocks
							for(int ii=0;ii<listBlocks.size()-1;ii++) {
								// AlignmentBlock ab=listBlocks.get(ii);

								AlignmentBlock ab1=listBlocks.get(ii);
								AlignmentBlock ab2=listBlocks.get(ii+1);

								int blockFrom=ab1.getReferenceStart();
								int blockTo=ab1.getReferenceStart()+ab1.getLength();
								int blockShouldBeInChunk=blockFrom/chunkSize;

								int inbetweenFrom=blockTo + 1;  // + 1; Assuming inclusive coordinates
								int inbetweenTo=ab2.getReferenceStart() - 1;  // -1; Assuming inclusive coordinates
								int inbetweenShouldBeInChunk=inbetweenFrom/chunkSize;

								// // CHECK: This has been commented out since before changes.
								// // 		 Do we need it for anything? //AB
								//If this is the first read we see, start chunking from here
								/*if(currentSeq.equals("")) {
									currentSeq=blockSource;
									saveAndSetCurrentChunk(blockSource, shouldBeInChunk);
								} else*/ 

								// Add alignment block to leftovers if not in right
								// chunk, else add to current chunk
								if((!currentSeq.equals(blockSource) || 
								   currentChunk!=blockShouldBeInChunk)) {
									addRegion(mapCellLeftoverAlignmentBlocks, barcodeIndex, 
												blockFrom, blockTo);
								} else {
									addRegion(mapCellAlignmentBlocks, barcodeIndex, 
												blockFrom, blockTo);
								}

								// Add inbetween block to leftovers if not in right
								// chunk, else add to current chunk
								if((!currentSeq.equals(blockSource) || 
								   currentChunk!=inbetweenShouldBeInChunk)) {
									addRegion(mapCellLeftoverInbetweens, barcodeIndex, 
												inbetweenFrom, inbetweenTo);
								} else {
									addRegion(mapCellInbetweens, barcodeIndex, 
												inbetweenFrom, inbetweenTo);
								}


							}
							// Last alignment block outside loop to stay in
							// bounds since inbetween does not exist here
							int ii = listBlocks.size()-1; // 0-based indexing so 
														  // last is size - 1 
							AlignmentBlock ab1=listBlocks.get(ii);

							int blockFrom=ab1.getReferenceStart();
							int blockTo=ab1.getReferenceStart()+ab1.getLength();
							int blockShouldBeInChunk=blockFrom/chunkSize;

							// // CHECK: This has been commented out since before changes.
							// // 		 Do we need it for anything? //AB
							//If this is the first read we see, start chunking from here
							/*if(currentSeq.equals("")) {
								currentSeq=blockSource;
								saveAndSetCurrentChunk(blockSource, shouldBeInChunk);
							} else*/ 

							// Add alignment block if in correct chunk
							if((!currentSeq.equals(blockSource) || 
							   currentChunk!=blockShouldBeInChunk)) {
								addRegion(mapCellLeftoverAlignmentBlocks, barcodeIndex, 
											blockFrom, blockTo);
							} else {
								addRegion(mapCellAlignmentBlocks, barcodeIndex, 
											blockFrom, blockTo);
							}

							// Here is more appropriate than after each alignment
							// block in read, since they are part of same
							// SAM record, which is what the counter
							// readRecords is keeping track of //AB
							keptRecords++;
						}
					} else {
						skippedBadUMI++;
					}
				} else {
					skippedWrongBC++;
				}
			}


		// BAM file is bulk type
		// (bamType.equals("--bulk")) must be, since only two options
		} else {  
	
			//Loop through all SAM records
			for (final SAMRecord samRecord : reader) {
				
				//Update user about progress
				readRecords++;
				if(readRecords%1000000 == 0){
					//Calculate progress
					int passedChunks=0;
					for(String seq:mapChunkStarts.keySet()) {
						if(seq.equals(currentSeq))
							break;
						passedChunks+=mapChunkStarts.get(seq).length;
					}
					passedChunks+=currentChunk;
					int prc=100*passedChunks/totalChunks;
					
					LogUtil.formatColumns(System.out, 25,
							prc+"%",
							"Kept/Read: "+keptRecords+"/"+readRecords,
							"@sequence: "+currentSeq,
							"WrongBC: "+skippedWrongBC,
							"badUMI: "+skippedBadUMI,
							"SkipDup: "+skippedDup,
							"noAlignBlocks: "+skippedNoAlignmentBlocks);
				}
					
				// Bulk verion; No check against BCs and UMIs

				// //Get UMI and BC for this read
				// String bcCellCurrentUMI=(String)samRecord.getAttribute("UB");
				// String bcCellCurrentCellBarcode=(String)samRecord.getAttribute("CB");
				
				// //Check if this is a cell to count
				// if(mapBarcodeIndex.keySet().contains(bcCellCurrentCellBarcode)) {
				// 	//If the read has no UMI nor BC then ignore it
				// 	if(bcCellCurrentCellBarcode!=null) {
				// 		//Check if duplicate read, if UMI present; ATAC, dedup by coordinate?
				// 		if(bcCellCurrentUMI!=null && 
				// 			bcCellCurrentUMI.equals(bcCellPreviousUMI) && 
				// 			bcCellCurrentCellBarcode.equals(bcCellPreviousCB)) {
				// 			skippedDup++;
				// 		} else {
				// 			//Remember for later
				// 			bcCellPreviousUMI=bcCellCurrentUMI;
				// 			bcCellPreviousCB=bcCellCurrentCellBarcode;



							// // Which cell is this?
							// int barcodeIndex=mapBarcodeIndex.get(bcCellCurrentCellBarcode);
							//
							// Bulk version; Only one barcode
							int barcodeIndex=0;


							//A read may have been split into multiple blocks.
							List<AlignmentBlock> listBlocks=samRecord.getAlignmentBlocks();
							if (listBlocks.size() == 0) { 
								skippedNoAlignmentBlocks++;
								continue;  	// Wouldn't it be a good idea to use
										   	// more continue statements instead
											// of all the current nested if else?
											// if .. continue combo doesn't
											// introduce deeper nesting.
											// //AB
							}

							String blockSource=samRecord.getContig();

							// If read starts outside chunk, change chunk,
							// since BAM is position sorted on start position
							// of reads.
							// "leftmost coordinates"
							// http://www.htslib.org/doc/samtools-sort.html
							int aa=listBlocks.get(0).getReferenceStart();
							int shouldBeInChunk=aa/chunkSize;
							if((!currentSeq.equals(blockSource) || 
							   currentChunk!=shouldBeInChunk)) {
								saveAndSetCurrentChunk(blockSource, shouldBeInChunk);
							}

							// Add the alignment blocks
							for(int ii=0;ii<listBlocks.size()-1;ii++) {
								// AlignmentBlock ab=listBlocks.get(ii);

								AlignmentBlock ab1=listBlocks.get(ii);
								AlignmentBlock ab2=listBlocks.get(ii+1);

								int blockFrom=ab1.getReferenceStart();
								int blockTo=ab1.getReferenceStart()+ab1.getLength();
								int blockShouldBeInChunk=blockFrom/chunkSize;

								int inbetweenFrom=blockTo + 1;  // + 1; Assuming inclusive coordinates
								int inbetweenTo=ab2.getReferenceStart() - 1;  // -1; Assuming inclusive coordinates
								int inbetweenShouldBeInChunk=inbetweenFrom/chunkSize;

								// // CHECK: This has been commented out since before changes.
								// // 		 Do we need it for anything? 
								//If this is the first read we see, start chunking from here
								/*if(currentSeq.equals("")) {
									currentSeq=blockSource;
									saveAndSetCurrentChunk(blockSource, shouldBeInChunk);
								} else*/ 

								// Add alignment block to leftovers if not in right
								// chunk, else add to current chunk
								if((!currentSeq.equals(blockSource) || 
								   currentChunk!=blockShouldBeInChunk)) {
									addRegion(mapCellLeftoverAlignmentBlocks, barcodeIndex, 
												blockFrom, blockTo);
								} else {
									addRegion(mapCellAlignmentBlocks, barcodeIndex, 
												blockFrom, blockTo);
								}

								// Add inbetween block to leftovers if not in right
								// chunk, else add to current chunk
								if((!currentSeq.equals(blockSource) || 
								   currentChunk!=inbetweenShouldBeInChunk)) {
									addRegion(mapCellLeftoverInbetweens, barcodeIndex, 
												inbetweenFrom, inbetweenTo);
								} else {
									addRegion(mapCellInbetweens, barcodeIndex, 
												inbetweenFrom, inbetweenTo);
								}


							}
							// Last alignment block outside loop to stay in 
							// bounds since inbetween does not exist here
							int ii = listBlocks.size()-1; // 0-based indexing so 
														  // last is size - 1 
							AlignmentBlock ab1=listBlocks.get(ii);

							int blockFrom=ab1.getReferenceStart();
							int blockTo=ab1.getReferenceStart()+ab1.getLength();
							int blockShouldBeInChunk=blockFrom/chunkSize;

							// // CHECK: This has been commented out since before changes.
							// // 		 Do we need it for anything? 
							//If this is the first read we see, start chunking from here
							/*if(currentSeq.equals("")) {
								currentSeq=blockSource;
								saveAndSetCurrentChunk(blockSource, shouldBeInChunk);
							} else*/ 

							// Add alignment block if in correct chunk
							if((!currentSeq.equals(blockSource) || 
							   currentChunk!=blockShouldBeInChunk)) {
								addRegion(mapCellLeftoverAlignmentBlocks, barcodeIndex, 
											blockFrom, blockTo);
							} else {
								addRegion(mapCellAlignmentBlocks, barcodeIndex, 
											blockFrom, blockTo);
							}

							// Here is more appropriate than after each alignment
							// block in read, since they are part of same
							// SAM record, which is what the counter
							// readRecords is keeping track of //AB
							keptRecords++;
				//	    }
				// 	} else {
				// 		skippedBadUMI++;
				// 	}
				// } else {
				// 	skippedWrongBC++;
				// }
			}

		}  // End of bulk file else



		//Get out the last ones from memory. TODO remember to fix this in countBAM too
		//
		// Is this TODO still relevant? 
		// Was here before I started changing the code. //AB
		saveCurrentChunk();
		reader.close();
		
		LogUtil.formatColumns(System.out, 25,
				"Done",
				"Kept/Read: "+keptRecords+"/"+readRecords,
				"@sequence: "+currentSeq,
				"WrongBC: "+skippedWrongBC,
				"badUMI: "+skippedBadUMI,
				"SkipDup: "+skippedDup);
	}
	
	
	
	/**
	 * Close the file
	 */
	public void close() throws IOException {
		if(raf!=null) {
			raf.close();
			raf=null;
		}
	}
	
	
	
	/////////////////////////////////////////////////////////////////////////////////////
	///////////////////// Read CellPiles    /////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////


	/**
	 * Read a string from file
	 */
	public String readString() throws IOException {
		int len=raf.readShort();
		StringBuilder sb=new StringBuilder();
		for(int i=0;i<len;i++) {
			sb.append(raf.readChar());
		}
		return sb.toString();
	}
	
	
	/**
	 * Open a file for reading
	 */
	public static CellPileFile open(File fCellpile) throws IOException {
		CellPileFile p=new CellPileFile();
		p.raf=new RandomAccessFile(fCellpile, "r");
		p.readHeader();
		return p;
	}

	/**
	 * Read the header of a file
	 */
	private void readHeader() throws IOException {
		//Go to the beginning of the file
		raf.seek(0);
		
		//Check that this is the correct type of file
		char[] compTo="CELLPILE".toCharArray();
		byte[] magicString=new byte["CELLPILE".length()]; //CELLPILE
		raf.readFully(magicString);
		for(int i=0;i<compTo.length;i++) {
			if(compTo[i]!=magicString[i]) {
				throw new RuntimeException("This is not a cellpile file");
			}
		}
		
		//Read version
		byte version=raf.readByte();
		System.out.println("Cellpile version: "+version);
		if(version!=2) {
			throw new RuntimeException("Unsupported version");
		}
		
		chunkSize=raf.readInt();
		
		//Read barcodes
		int numBarcodes=raf.readInt();
		listBarcodes=new ArrayList<String>();
		for(int i=0;i<numBarcodes;i++) {
			listBarcodes.add(readString());
		}
		prepareBarcodeIndexLookup();

		//Read the names and positions of all chunks
		int numSequences=raf.readInt();
		//System.out.println("#seq: "+numSequences);
		for(int i=0;i<numSequences;i++) {
			String nameOfChunk=readString();
			int numChunk=raf.readInt();
			long[] chunkpos=new long[numChunk];
			for(int curc=0;curc<chunkpos.length;curc++) {
				chunkpos[curc]=raf.readLong();
			}
			mapChunkStarts.put(nameOfChunk, chunkpos);
		}
	}
	
	
	/**
	 * Convert clusters naming cells by string, to IDs  (list of lists)
	 */
	public int[][] convertBarcodeNamesToIDs(String[][] bclistS){
		int[][] outlist=new int[bclistS.length][];
		for(int i=0;i<bclistS.length;i++) {
			outlist[i]=convertBarcodeNamesToIDs(bclistS[i]);
		}
		return outlist;
	}
	
	
	/**
	 * Convert clusters naming cells by string, to IDs
	 */
	public int[] convertBarcodeNamesToIDs(String[] bclistS){
		int[] list2=new int[bclistS.length];
		for(int j=0;j<bclistS.length;j++) {
			list2[j] = mapBarcodeIndex.get(bclistS[j]);
		}
		return list2;
	}
	
	
	/**
	 * Build pileups for each given cell group
	 * @param windowSeq     Sequence to look up
	 * @param windowFrom    From-position
	 * @param windowTo      To-position
	 * @param numdiv        Number of subdivisions. Minimum 2
	 * @param cellGroups    List of lists of barcodes (the clusters to pile-up for)
	 */
	public Pileup buildPileup(
			String windowSeq, int windowFrom, int windowTo, int numdiv,
			int[][] cellGroups, String[] clusterNames) throws IOException {

		//Create the basic return object
		Pileup pileup=new Pileup();
		pileup.seq=windowSeq;
		pileup.from=windowFrom;
		pileup.to=windowTo;
		pileup.numdiv=numdiv;
		pileup.cellCluster=cellGroups;
		pileup.clusterNames=clusterNames;

		//Set the cell counts
		pileup.clusterCellCount=new int[clusterNames.length];
		for(int i=0;i<cellGroups.length;i++)
			pileup.clusterCellCount[i] = cellGroups[i].length;
			
		//Pileup grid spacing
		double dx=(windowTo-windowFrom+1)/(double)(numdiv-1);
		
		//Do some checking
		if(numdiv<2)
			throw new RuntimeException("Too few subdivisions");
		if(mapChunkStarts.get(windowSeq)==null) {
			throw new RuntimeException("Sequence does not exist: "+windowSeq);
			
		}

		//Quick lookup, which cell on which track
		TreeMap<Integer, Integer> mapBarcodeTrack=new TreeMap<Integer, Integer>();
		for(int curTrack=0;curTrack<cellGroups.length;curTrack++)
			for(int bc:cellGroups[curTrack])
				mapBarcodeTrack.put(bc,curTrack);
		
		//Allocate the pileup tracks
		int[][] outTracksAlignmentBlocks=new int[cellGroups.length][numdiv];
		int[][] outTracksInbetweens=new int[cellGroups.length][numdiv];

		//Figure out chunks to read
		int chunkFrom=windowFrom/chunkSize;
		int chunkTo=windowTo/chunkSize+1;
		chunkFrom = PileUtil.clamp(chunkFrom, 0, mapChunkStarts.get(windowSeq).length-1);
		chunkTo = PileUtil.clamp(chunkTo, 0, mapChunkStarts.get(windowSeq).length-1);
		// Always read one more chunk upstream of range user asks for.
		// This guarantees that alignmentblocks starting in upstream chunk
		// and ending in current chunk are included.
		// Exception is if chunkFrom is 0, then there is no upstream chunk.
		if (chunkFrom != 0) {
			chunkFrom = chunkFrom - 1;
		}

		//Iterate through all the chunks
		int counted=0;
		for(int curChunk=chunkFrom;curChunk<=chunkTo;curChunk++) {

			long chunkPos=mapChunkStarts.get(windowSeq)[curChunk];
			if(chunkPos!=0) {

				raf.seek(chunkPos);

				readRegionType(outTracksAlignmentBlocks, windowFrom, windowTo,
								numdiv, counted, mapBarcodeTrack, dx);
				
				readRegionType(outTracksInbetweens, windowFrom, windowTo,
								numdiv, counted, mapBarcodeTrack, dx);
			}
		}
		
		pileup.alignmentBlockTracks=outTracksAlignmentBlocks;
		pileup.inbetweenTracks=outTracksInbetweens;

		return pileup;
	}
	

	// Reads one region type (alignment blocks or inbetweens) from cellpilefile.
	// The region types are always stacked after each other in each chunk,
	// so can call this twice in a row to get both.
	public void readRegionType(int[][] tracks, int windowFrom, int windowTo,
								int numdiv, int counted,
								TreeMap<Integer, Integer> mapBarcodeTrack,
								double dx) throws IOException {

		//Loop through all cells represented in this chunk
		int numCells=raf.readInt();

		nextcell: for(int curCell=0;curCell<numCells;curCell++) {
			//Filter cells, or select the right track
			int cellID=raf.readInt();   

			Integer toTrack=mapBarcodeTrack.get(cellID);
			if(toTrack!=null) {
				//Check how many regions
				int[] thisTrack=tracks[toTrack];
				int numRegions=raf.readShort();
				
				//Read each region
				for(int curRegion=0;curRegion<numRegions;curRegion++) {
					int posLeft=raf.readInt();
					int posRight=raf.readInt();

					//We may be able to quit early if lucky
					if(posLeft<=windowTo) {
						//Transform to pileup coordinates
						posLeft=(int)((posLeft-windowFrom)/dx);
						posRight=(int)((posRight-windowFrom)/dx);
						
						//Restrict to only count within track limits
						posLeft=Math.max(posLeft,0);
						posRight=Math.min(posRight, numdiv-1);
						
						for(int j=posLeft;j<=posRight;j++) {
							thisTrack[j]++;
							counted++;
						}
					} else {
						//Get this to work after basic edition working - dangerous!
						//raf.skipBytes(4*(numRegions-1-i));
						//continue nextcell;
					}
				}
			}
		}
	}



	/**
	 * Get a list of sequences represented
	 */
	public Collection<String> getListSequences(){
		return Collections.unmodifiableSet(mapChunkStarts.keySet());
	}
	public String[] getListSequencesAsArray(){
		return new ArrayList<String>(mapChunkStarts.keySet()).toArray(new String[0]);
	}

	/**
	 * Get a list of barcodes represented
	 */
	public List<String> getListBarcodes(){
		return Collections.unmodifiableList(listBarcodes);
	}
	public String[] getListBarcodesAsArray(){
		return listBarcodes.toArray(new String[0]);
	}

	public static void printInfo(File fCellpile) throws IOException {

		RandomAccessFile raf=new RandomAccessFile(fCellpile, "r");

		//Go to the beginning of the file
		raf.seek(0);

		//Check that this is the correct type of file
		char[] compTo="CELLPILE".toCharArray();
		byte[] magicString=new byte["CELLPILE".length()]; //CELLPILE
		raf.readFully(magicString);
		for(int i=0;i<compTo.length;i++) {
			if(compTo[i]!=magicString[i]) {
				throw new RuntimeException("This is not a cellpile file");
			}
		}

		//Read version
		byte version=raf.readByte();
		System.out.println("Cellpile version: "+version);

	}
}
