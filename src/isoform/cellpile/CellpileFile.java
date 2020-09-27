package isoform.cellpile;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
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
 * Generate a .cellpile from a BAM-file
 * 
 * 
 * Consider swapping to a buffered implementation at least for the read interface
 * 
 * @author Johan Henriksson
 *
 */
public class CellpileFile {
	
	private RandomAccessFile raf;
	private int chunkSize=10000; //10kb


	private TreeMap<String, long[]> mapChunkStarts=new TreeMap<String, long[]>();
	private HashMap<String, Integer> mapBarcodeIndex=new HashMap<String, Integer>();
	private ArrayList<String> listBarcodes;
	
	private int dups=0;
	private int kept=0;
	
	public ArrayList<String> getListBarcodes(){
		return new ArrayList<String>(listBarcodes);
	}

	/**
	 * Write a string to file
	 */
	public void writeString(String s) throws IOException {
		raf.writeShort((short)s.length());
		for(char c:s.toCharArray())
			raf.writeChar(c);
	}
	
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
	 * Constructor - cannot be called from the outside
	 */
	private CellpileFile() {
	}
	

	/**
	 * Read chromosome sizes and calculate chunk sizes
	 */
	private void readChromosomeSizes(File f) throws IOException {
		
		
		//samtools faidx input.fa
		//cut -f1,2 input.fa.fai > sizes.genome
		
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
	public static CellpileFile writeFile(
			File fCellpile, 
			File fChromSizes,
			File fBAM,
			ArrayList<String> listBarcodes) throws IOException {
		
		CellpileFile p=new CellpileFile();
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

		//Calculate where this is and store it for later
		//long endOfHeader=raf.getFilePointer();
		
	}
	
	
	//Current set of regions that we are collecting
	private TreeMap<Integer, IntArrayList> mapCellRegions=new TreeMap<Integer, IntArrayList>();
	private String currentSeq;
	private int currentChunk;

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
	}
	

	/**
	 * Save chunk... if there is anything to save
	 */
	private void saveCurrentChunk() throws IOException {

		if(!mapCellRegions.isEmpty()) {
			//Store pointer to where this chunk starts in the file
			mapChunkStarts.get(currentSeq)[currentChunk]=raf.getFilePointer();

			//Store the chunk
			raf.writeInt(mapCellRegions.size());
			for(int bc:mapCellRegions.keySet()) {
				raf.writeInt(bc);
				IntArrayList ia=mapCellRegions.get(bc);
				int[] list=ia.elements(); //Note, no copying, and wrong size
				int len=ia.size();
				for(int i=0;i<len;i++) {
					raf.writeInt(list[i]);
				}
			}
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


		
		System.out.println("Moved to chunk "+newSeq+"\t"+newChunkNum);
	}
	
	
	/**
	 * Perform the counting from a BAM file. Can be called multiple times
	 */
	public void countReads(File fBAM) throws IOException {
		final SamReader reader = SamReaderFactory.makeDefault().open(fBAM);
		
		//Set initial chunk position to nowhere
		currentSeq="";
		currentChunk=0;

		//To know progress,
		int readRecords=0;

		//Store status for de-duplication
		SAMRecord previousRecord=null;
		String bcCellPreviousU=null;
		String bcCellPreviousC=null;
				
		for (final SAMRecord samRecord : reader) {
			readRecords++;
			if(readRecords%10000 == 0){
				System.out.println("Kept/Read: "+kept+"/"+readRecords);
			}
			/*if(readRecords%100000 == 0){
				break;
			}*/
				
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
								
								int shouldBeInChunk=blockFrom/chunkSize;
								
								//If this is the first read we see, start chunking from here
								if(currentSeq.equals("")) {
									currentSeq=blockSource;
									saveAndSetCurrentChunk(blockSource, shouldBeInChunk);
								} else if(curAB==0 && (!currentSeq.equals(blockSource) || currentChunk!=shouldBeInChunk)) {
									//Check if we are still within the same chunk. otherwise move on.
									saveAndSetCurrentChunk(blockSource, shouldBeInChunk);
								}
								
								addRegion(barcodeIndex, blockFrom, blockTo);
								kept++;
							}
						}
					}
				}
			} else {
				previousRecord=samRecord;
			}
		}
		
		//Get out the last ones from memory. TODO remember to fix this in countBAM too
		saveCurrentChunk();
		
		System.out.println("Kept: "+kept);
		System.out.println("Duplicates: "+dups);

		reader.close();
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
	
	/**
	 * Finalizer. Ensures somewhat fast closing of the file - although this is better done explicitly!
	 */
	protected void finalize() throws Throwable {
		close();
	}

	
	
	/////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////


	/**
	 * Open a file for reading
	 */
	public static CellpileFile open(File fCellpile) throws IOException {
		CellpileFile p=new CellpileFile();
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
		System.out.println("#seq: "+numSequences);
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
	 * Convert clusters naming cells by string, to IDs
	 */
	public int[][] convertBarcodeNamesToIDs(String[][] bclistS){
		int[][] outlist=new int[bclistS.length][];
		for(int i=0;i<bclistS.length;i++) {
			String[] list1=bclistS[i];
			int[] list2=new int[list1.length];
			for(int j=0;j<list1.length;j++) {
				list2[j] = mapBarcodeIndex.get(list1[j]);
			}
			outlist[i]=list2;
		}
		return outlist;
	}
	
	
	/**
	 * Build pileups for each given cell group
	 * @param seq     Sequence to look up
	 * @param from    From-position
	 * @param to      To-position
	 * @param numdiv  Number of subdivisions. Minimum 2
	 * @param cellGroups  List of lists of barcodes (the clusters to pile-up for)
	 */
	public int[][] buildPileup(
			String seq, int from, int to, 
			int numdiv,
			int[][] cellGroups) throws IOException {
		
		//Pileup grid spacing
		double dx=(to-from+1)/(double)(numdiv-1);
		
		//Do some checking
		if(numdiv<2)
			throw new RuntimeException("Too few subdivisions");
		if(mapChunkStarts.get(seq)==null) {
			System.out.println(mapChunkStarts.keySet());
			throw new RuntimeException("Sequence does not exist: "+seq);
			
		}

		//Quick lookup, which cell on which track
		TreeMap<Integer, Integer> mapBarcodeTrack=new TreeMap<Integer, Integer>();
		for(int curTrack=0;curTrack<cellGroups.length;curTrack++)
			for(int bc:cellGroups[curTrack])
				mapBarcodeTrack.put(bc,curTrack);
		
		//Allocate the pileup tracks
		int[][] outTracks=new int[cellGroups.length][numdiv];

		//Iterate through all the chunks
		int chunkFrom=from/chunkSize;
		int chunkTo=to/chunkSize+1;
		for(int curChunk=chunkFrom;curChunk<=chunkTo;curChunk++) {
			long chunkPos=mapChunkStarts.get(seq)[curChunk];
			if(chunkPos!=0) {
				raf.seek(chunkPos);
				
				//Loop through all cells represented in this chunk
				int numCells=raf.readInt();
				nextcell: for(int curCell=0;curCell<numCells;curCell++) {
					//Filter cells, or select the right track
					int cellID=raf.readInt();
					Integer toTrack=mapBarcodeTrack.get(cellID);
					if(toTrack!=null) {
						int[] thisTrack=outTracks[toTrack];
						int numRegions=raf.readShort();
						
						for(int curRegion=0;curRegion<numRegions;curRegion++) {
							int posLeft=raf.readInt();
							int posRight=raf.readInt();

							//We may be able to quit early if lucky
							if(posLeft<=to) {
								//Transform to pileup coordinates
								posLeft=(int)((posLeft-from)/dx);
								posRight=(int)((posRight-from)/dx);
								
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
		
		return outTracks;
	}
	
	
	
}
