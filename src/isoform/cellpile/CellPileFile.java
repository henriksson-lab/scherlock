package isoform.cellpile;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.StringTokenizer;
import java.util.TreeMap;

import cern.colt.list.IntArrayList;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * Interface for CellPile files. Able to generate from BAM-files
 * 
 * Consider swapping to a buffered implementation for the read interface (need to verify if this improves speed)
 * 
 * @author Johan Henriksson
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
	///////////////////// Reading CellPiles /////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////

	/////////////////////
	//For building: Current set of regions that we are collecting
	private TreeMap<Integer, IntArrayList> mapCellRegions=new TreeMap<Integer, IntArrayList>();
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

		p.countReads(fBAM);
		
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
		byte version=1;
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
	private void addRegion(int bcid, int from, int to) {
		IntArrayList regs=mapCellRegions.get(bcid);
		if(regs==null) {
			mapCellRegions.put(bcid, regs=new IntArrayList(1000));
		}
		regs.add(from);
		regs.add(to);
//		System.out.println("store "+bcid+"   "+listBarcodes.get(bcid)+"   "+from+"    "+to);
	}
	

	/**
	 * Save current chunk
	 */
	private void saveCurrentChunk() throws IOException {
		//Only store chunks with data. Pointer is default 0=missing
		if(!mapCellRegions.isEmpty()) {
			//Ensure that this chunk has been allocated in the table; otherwise ignore it
			long[] ptrs=mapChunkStarts.get(currentSeq);
			if(ptrs!=null) {
				if(currentChunk<ptrs.length) {
					//Store pointer to where this chunk starts in the file
					ptrs[currentChunk]=raf.getFilePointer();
//					System.out.println("file pointer "+raf.getFilePointer());

					//Store the chunk
					int numCells=mapCellRegions.size();
//					System.out.println("Storing # cells "+numCells);
					raf.writeInt(numCells);
					for(int cellBarcode:mapCellRegions.keySet()) {
						raf.writeInt(cellBarcode);
						//System.out.println("---cbi "+cellBarcode);
						IntArrayList ia=mapCellRegions.get(cellBarcode);
						int[] list=ia.elements(); //Note, this is not a copy of the array, and it is the wrong size
						int numRegion2=ia.size();
						raf.writeShort(numRegion2/2);
						for(int i=0;i<numRegion2;i++) {
							raf.writeInt(list[i]);
						}
					}
					
					//System.out.println("Stored chunk, "+currentChunk+ " #cells "+mapCellRegions.size());
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
	
	/**
	 * Save current chunk and move to the next position
	 */
	private void saveAndSetCurrentChunk(String newSeq, int newChunkNum) throws IOException {
		saveCurrentChunk();
		
		//Reset chunk store so we can fill with the next region
		mapCellRegions.clear();
		currentSeq=newSeq;
		currentChunk=newChunkNum;
	}
	
	

	/**
	 * Perform the counting from a BAM file. Can be called multiple times
	 */
	public void countReads(File fBAM) throws IOException {
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

		//This is to keep track of duplicates.
		//Approximate, as the same cDNA can be fragmented multiple times in the library prep
		String bcCellPreviousUMI="";
		String bcCellPreviousCB="";
				
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
				System.out.println(""+prc+"%\tKept/Read: "+keptRecords+"/"+readRecords+"\tOn sequence: +"+currentSeq);
			}
				
			//Get UMI and BC for this read
			String bcCellCurrentUMI=(String)samRecord.getAttribute("UB");
			String bcCellCurrentCellBarcode=(String)samRecord.getAttribute("CB");
			
			//Check if this is a cell to count
			if(mapBarcodeIndex.keySet().contains(bcCellCurrentCellBarcode)) {
				//If the read has no UMI nor BC then ignore it
				if(bcCellCurrentUMI!=null && bcCellCurrentCellBarcode!=null) {
					//Check if duplicate read
					if(bcCellCurrentUMI.equals(bcCellPreviousUMI) && bcCellCurrentCellBarcode.equals(bcCellPreviousCB)) {
						//Do nothing, ignore read
						//System.out.println("Got duplicate");
					} else {
						//Remember for later
						bcCellPreviousUMI=bcCellCurrentUMI;
						bcCellPreviousCB=bcCellCurrentCellBarcode;

						//Which cell is this?
						int barcodeIndex=mapBarcodeIndex.get(bcCellCurrentCellBarcode);
						
						//A read may have been split into multiple blocks. 
						//Count these separately. Naive assumption that these are split over introns... is this correct?
						List<AlignmentBlock> listBlocks=samRecord.getAlignmentBlocks();
						//System.out.println("#alignment blocks "+listBlocks.size());
						for(int curAB=0;curAB<listBlocks.size();curAB++) {
							AlignmentBlock ab=listBlocks.get(curAB);
															
							String blockSource=samRecord.getContig();
							int blockFrom=ab.getReferenceStart();
							int blockTo=ab.getReferenceStart()+ab.getLength();
							
							int shouldBeInChunk=blockFrom/chunkSize;
							
							//If this is the first read we see, start chunking from here
							/*if(currentSeq.equals("")) {
								currentSeq=blockSource;
								saveAndSetCurrentChunk(blockSource, shouldBeInChunk);
							} else*/ 
							if(curAB==0 && (!currentSeq.equals(blockSource) || currentChunk!=shouldBeInChunk)) {
								//Check if we are still within the same chunk. otherwise move on.
								saveAndSetCurrentChunk(blockSource, shouldBeInChunk);
							}
							
							addRegion(barcodeIndex, blockFrom, blockTo);
							keptRecords++;
						}
					}
				} else {
					//System.out.println("Incomplete BAM record");
				}
			} else {
//				System.out.println("Skipping cell: "+bcCellCurrentCellBarcode);
			}
		}
		
		//Get out the last ones from memory. TODO remember to fix this in countBAM too
		saveCurrentChunk();
		reader.close();
		
		System.out.println("Final - Kept/Read: "+keptRecords+"/"+readRecords);
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
	///////////////////// Writing CellPiles /////////////////////////////////////////////
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
		System.out.println("Version: "+version);
		if(version!=1) {
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
	//			if(chunkpos[curc]!=0)
//					System.out.println("got a chunk");
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
			String windowSeq, int windowFrom, int windowTo, 
			int numdiv,
			int[][] cellGroups, String[] clusterNames) throws IOException {
		
		//Create the basic return object
		Pileup pileup=new Pileup();
		pileup.seq=windowSeq;
		pileup.from=windowFrom;
		pileup.to=windowTo;
		pileup.numdiv=numdiv;
		pileup.cellGroups=cellGroups;
		pileup.clusterNames=clusterNames;
		
		
		//Pileup grid spacing
		double dx=(windowTo-windowFrom+1)/(double)(numdiv-1);
		
		//Do some checking
		if(numdiv<2)
			throw new RuntimeException("Too few subdivisions");
		if(mapChunkStarts.get(windowSeq)==null) {
			//System.out.println(mapChunkStarts.keySet());
			throw new RuntimeException("Sequence does not exist: "+windowSeq);
			
		}

		//Quick lookup, which cell on which track
		TreeMap<Integer, Integer> mapBarcodeTrack=new TreeMap<Integer, Integer>();
		for(int curTrack=0;curTrack<cellGroups.length;curTrack++)
			for(int bc:cellGroups[curTrack])
				mapBarcodeTrack.put(bc,curTrack);
		
		//Allocate the pileup tracks
		int[][] outTracks=new int[cellGroups.length][numdiv];

		//Iterate through all the chunks
		int chunkFrom=windowFrom/chunkSize;
		int chunkTo=windowTo/chunkSize+1;
		for(int curChunk=chunkFrom;curChunk<=chunkTo;curChunk++) {
			long chunkPos=mapChunkStarts.get(windowSeq)[curChunk];
			if(chunkPos!=0) {
				raf.seek(chunkPos);
//				System.out.println("go to fp "+chunkPos);
				
				//Loop through all cells represented in this chunk
				int numCells=raf.readInt();
	//			System.out.println("check chunk "+curChunk+"  "+numCells);
				nextcell: for(int curCell=0;curCell<numCells;curCell++) {
					//Filter cells, or select the right track
					int cellID=raf.readInt();                                    ///////////// value totally in the blue
					Integer toTrack=mapBarcodeTrack.get(cellID);
					if(toTrack!=null) {
//						System.out.println("---cbi "+cellID+"   "+toTrack);
						//Check how many regions
			//			System.out.println("cell "+cellID);
						int[] thisTrack=outTracks[toTrack];
						int numRegions=raf.readShort();
						
						//Read each region
						for(int curRegion=0;curRegion<numRegions;curRegion++) {
							int posLeft=raf.readInt();
							int posRight=raf.readInt();


							//We may be able to quit early if lucky
							if(posLeft<=windowTo) {
								//System.out.println("---cbi "+cellID+"   "+toTrack+"   "+curRegion+"\t"+posLeft+"\t"+posRight);
								//Transform to pileup coordinates
								posLeft=(int)((posLeft-windowFrom)/dx);
								posRight=(int)((posRight-windowFrom)/dx);
								
								//Restrict to only count within track limits
								posLeft=Math.max(posLeft,0);
								posRight=Math.min(posRight, numdiv-1);
								
								for(int j=posLeft;j<=posRight;j++) {
									thisTrack[j]++;
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
		}
		
		pileup.tracks=outTracks;
		return pileup;
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
	
	
	
	
}
