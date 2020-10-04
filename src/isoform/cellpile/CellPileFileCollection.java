package isoform.cellpile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.TreeSet;

import cern.colt.list.IntArrayList;

/**
 * A meta interface to work with several piles in parallel. This allows 
 * 
 * @author Johan Henriksson
 *
 */
public class CellPileFileCollection {

	public ArrayList<CellPileFile> listFile=new ArrayList<CellPileFile>();
	public ArrayList<String> listPileName=new ArrayList<String>();
	
	public void addFile(CellPileFile pile, String name, String bcPrefix) {
		listFile.add(pile);
		listPileName.add(name);
	}
	

	/**
	 * Convert barcodes to IDs for a given fileID
	 */
	private int[] convertBarcodeNamesToIDs(int forFile, String[] bclistS, String[] pilename, String thisPile){
		CellPileFile f=listFile.get(forFile);
		IntArrayList list2=new IntArrayList(bclistS.length);
		for(int j=0;j<bclistS.length;j++) {
			if(pilename[j].equals(thisPile)) {
				list2.add(f.mapBarcodeIndex.get(bclistS[j]));
			}
		}
		list2.trimToSize();
		return list2.elements();
	}
	
	
	/**
	 * Build a pileup that sums up all pileups
	 */
	public Pileup buildPileup(
			String windowSeq, int windowFrom, int windowTo, 
			int numdiv,
			String[] cellBC, String[] cellFile, String[] cellCluster) throws IOException {
		
		//Figure out what cluster names there are
		TreeSet<String> allClusterNames=new TreeSet<String>();
		for(String s:cellCluster)
			allClusterNames.add(s);
		
		//Assign numbers to each cluster
		TreeMap<String, Integer> mapClusterID=new TreeMap<String, Integer>();
		String[] listClusterNames=new String[allClusterNames.size()];
		int curAssId=0;
		for(String s:cellCluster) {
			mapClusterID.put(s,curAssId);
			listClusterNames[curAssId]=s;
			curAssId++;
		}
		int numClusters=curAssId;
		
		Pileup totalp=null;
		for(int i=0;i<listFile.size();i++) {
			
			//Generate barcode mappings for this file, which is a subset of all BCs given
			int[][] listBCs=new int[numClusters][];
			for(int j=0;j<numClusters;j++) {
				listBCs[j] = convertBarcodeNamesToIDs(j, cellBC, cellFile, listPileName.get(j));
			}
			
			//Perform the pilup for this file
			Pileup p=listFile.get(i).buildPileup(
					windowSeq, windowFrom, windowTo, 
					numdiv,
					listBCs, listClusterNames);
			
			if(i==0)
				totalp=p;
			else
				totalp.addPileup(p);
		}
		return totalp;
	}

	
}
