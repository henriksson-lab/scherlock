package isoform.cellpile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;

import cern.colt.list.IntArrayList;
import isoform.trackRenderer.Track;
import isoform.trackRenderer.TrackGTF;
import isoform.trackRenderer.TrackBed;
import isoform.trackRenderer.TrackPileup;
import isoform.trackRenderer.TrackRenderer;
import isoform.util.GtfParser;
import isoform.util.PileUtil;
import isoform.util.Range;

/**
 * High-level interface to work with several piles in parallel
 * 
 * @author Johan Henriksson and Anton Björk
 *
 */
public class CellPileManager {

	//// List of piles
	public ArrayList<CellPileFile> listPileFile=new ArrayList<CellPileFile>();
	public ArrayList<String> listPileName=new ArrayList<String>();
	
	//// List of tracks
	public ArrayList<Track> listTracks=new ArrayList<Track>();
	
	
	
	/**
	 * Add a cellpile file
	 */
	public void addCellPile(CellPileFile pile, String name) {
		listPileFile.add(pile);
		listPileName.add(name);
	}
	
	/**
	 * Get the number of cell piles
	 */
	public int getNumCellPiles() {
		return listPileFile.size();
	}

	/**
	 * Get the names of piles
	 */
	public List<String> getPileNames(){
		return listPileName;
	}
	
	/**
	 * Add a general track
	 */
	public void addTrack(String name, Track t) {
		t.trackName=name;
		listTracks.add(t);
	}
	
	/**
	 * Add a GTF track
	 */
	public void addTrackGTF(String name, File fname) throws IOException {
		GtfParser gtf=new GtfParser(fname);
		TrackGTF t=new TrackGTF(gtf);
		t.trackName=name;
		listTracks.add(t);
	}


	/**
	 * Parse and dump a gtf to python
	 */
	public GtfParser dumpTrackGTF(String name, File fname) throws IOException {
		GtfParser gtf=new GtfParser(fname);
		return gtf;
		// TrackGTF t=new TrackGTF(gtf);
		// t.trackName=name;
		// listTracks.add(t);
	}

	
	
	/**
	 * Add a Bed file track
	 */
	public void addTrackBed(String name, String fname) throws IOException {
		TrackBed t=new TrackBed(fname);
		t.trackName=name;
		listTracks.add(t);
	}

	
	
	
	/**
	 * Get the barcodes of one pileup
	 */
	public List<String> getBarcodes(int i){
		return listPileFile.get(i).getListBarcodes();
	}
	
	
	/**
	 * Convert barcodes to IDs for a given fileID
	 */
	private int[] convertBarcodeNamesToIDs(
			int forPile, 
			String[] cellBC, String[] cellPile, String[] cellCluster, 
			String forCluster){
		CellPileFile pile=listPileFile.get(forPile);
		if (pile == null) {
			throw new RuntimeException("No pile with id " + forPile);
		}
		String thisPileName=listPileName.get(forPile);
		
		
		IntArrayList list2=new IntArrayList(cellBC.length);
		for(int j=0;j<cellBC.length;j++) {
			if(cellPile[j].equals(thisPileName) && cellCluster[j].equals(forCluster)) {

				// Add to list if barcode is found, else crash because of mismatch in barcode lists
				// null if cellBC[j] not in pile.mapBarcodeIndex
				Integer aa = pile.mapBarcodeIndex.get(cellBC[j]);  
				if (aa != null) {  // Integer can be null but int cannot, so check before convert
					int bb = aa.intValue();
					list2.add(bb); // IntArrayList only takes int, not Integer
				} else {
// throw new RuntimeException("No barcode " + cellBC[j].toString() +	" in cellpile file. Aborting");
//TODO: better: store # missing barcodes in the output object. output this as statistics to the user
				}

				// Original line that was crashing with nullpointer exception sometimes
				// list2.add(pile.mapBarcodeIndex.get(cellBC[j]));
			}
		}
		list2.trimToSize();
		return list2.elements();
	}
	
	
	public static void print(String[] ss) {
		for(String s:ss)
			System.out.print(s+" ");
		System.out.println();
	}
	
	/**
	 * Build a pileup that sums up all pileups
	 */
	public Pileup buildPileup(
			String windowSeq, int windowFrom, int windowTo, 
			int numdiv,
			String[] cellBC, String[] cellFile, String[] cellCluster) throws IOException {
		
		
		if (cellBC == null) {
			throw new RuntimeException("cellBC is null!");
		}
		if (cellFile == null) {
			throw new RuntimeException("cellFile is null!");
		}
		if (cellCluster == null) {
			throw new RuntimeException("cellCluster is null!");
		}
		
		
		System.out.println("--- got BC");
		print(cellBC);
		System.out.println("--- got file");
		print(cellFile);
		System.out.println("--- got cluster");
		print(cellCluster);
		
		//Figure out what cluster names there are
		TreeSet<String> allClusterNames=new TreeSet<String>();
		for(String s:cellCluster)
			allClusterNames.add(s);
		int numClusters=allClusterNames.size();

		
		//Assign numbers to each cluster
		TreeMap<String, Integer> mapClusterID=new TreeMap<String, Integer>();
		String[] listClusterNames=new String[numClusters];
		int curAssId=0;
		for(String s:allClusterNames) {
			mapClusterID.put(s,curAssId);
			listClusterNames[curAssId]=s;
			curAssId++;
		}
		
		Pileup totalp=null;
		for(int iFile=0;iFile<listPileFile.size();iFile++) {
			
			//Generate barcode mappings for this file, which is a subset of all BCs given
			int[][] listBCs=new int[numClusters][];
			for(int iCluster=0;iCluster<numClusters;iCluster++) {
				listBCs[iCluster] = convertBarcodeNamesToIDs(iFile, cellBC, cellFile, 
					cellCluster, listClusterNames[iCluster]);
				System.out.println("got # cells: "+listBCs[iCluster].length);
			}
			System.out.println("------------------");
			
			//Perform the pilup for this file
			Pileup p=listPileFile.get(iFile).buildPileup(
					windowSeq, windowFrom, windowTo, 
					numdiv,
					listBCs, listClusterNames);
			
			if(iFile==0)
				totalp=p;
			else
				totalp.addPileup(p);
		}
		
		
		//Transfer all the tracks
		//for(TrackListener)
		
		return totalp;
	}

	
	/**
	 * Build a pileup, everything in one group
	 */
	public Pileup buildTotalPileup(
			String windowSeq, int windowFrom, int windowTo, 
			int numdiv) throws IOException {
		
		ArrayList<String> arrAllBC=new ArrayList<String>();
		ArrayList<String> arrAllFile=new ArrayList<String>();
		for(int i=0;i<listPileFile.size();i++) {
			CellPileFile cp=listPileFile.get(i);
			String pname=listPileName.get(i);
			
			List<String> list=cp.getListBarcodes();
			arrAllBC.addAll(list);
			for(int j=0;j<list.size();j++)
				arrAllFile.add(pname);
		}
		
		String[] cellBC=arrAllBC.toArray(new String[0]);
		String[] cellFile=arrAllFile.toArray(new String[0]);
		String[] cellCluster=PileUtil.getRepeatedString("", cellBC.length);
		
		return buildPileup(windowSeq, windowFrom, windowTo, numdiv, cellBC, cellFile, cellCluster);
	}
	
	
	/**
	 * Create a renderer. Pileup can be null
	 */
	public TrackRenderer render(String windowSeq, int windowFrom, int windowTo, 
				Pileup pileup, 
				boolean showInbetweens, boolean individualTrackScaling) {

		TrackRenderer renderer=new TrackRenderer(windowSeq, windowFrom, windowTo,
			individualTrackScaling);

		if(pileup!=null)
			renderer.addTrack(new TrackPileup(pileup, false, showInbetweens));
		for(Track t:listTracks)
			renderer.addTrack(t);
		
		return renderer;
	}
	
	

	/**
	 * Get the range of a gene. Gene can be symbol or ID
	 */
	public Range getRangeForGene(String gene) {
		
		
		for(Track t:listTracks) {
			if(t instanceof TrackGTF) {
				Range r=((TrackGTF)t).gtf.getRangeForGene(gene);
				if(r!=null)
					return r;
			}
		}
		return null;
	}

	
	
}
