package isoform.cellpile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
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
	
	public void addFile(CellPileFile pile, String name) {
		listFile.add(pile);
		listPileName.add(name);
	}
	
	public int getNumFiles() {
		return listFile.size();
	}

	/**
	 * Get the barcodes of one pileup
	 */
	public List<String> getBarcodes(int i){
		return listFile.get(i).getListBarcodes();
	}
	
	
	/**
	 * Convert barcodes to IDs for a given fileID
	 */
	private int[] convertBarcodeNamesToIDs(int forFile, String[] cellBC, String[] cellFile, String[] cellCluster, String thisCluster){
		CellPileFile f=listFile.get(forFile);
		String thisFileName=listPileName.get(forFile);
		IntArrayList list2=new IntArrayList(cellBC.length);
		for(int j=0;j<cellBC.length;j++) {
			if(cellFile[j].equals(thisFileName) && cellCluster[j].equals(thisCluster)) {
				list2.add(f.mapBarcodeIndex.get(cellBC[j]));
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
		for(int iFile=0;iFile<listFile.size();iFile++) {
			
			//Generate barcode mappings for this file, which is a subset of all BCs given
			int[][] listBCs=new int[numClusters][];
			for(int iCluster=0;iCluster<numClusters;iCluster++) {
				listBCs[iCluster] = convertBarcodeNamesToIDs(iFile, cellBC, cellFile, cellCluster, listClusterNames[iCluster]);
				System.out.println("got # cells: "+listBCs[iCluster].length);
			}
			System.out.println("------------------");
			
			//Perform the pilup for this file
			Pileup p=listFile.get(iFile).buildPileup(
					windowSeq, windowFrom, windowTo, 
					numdiv,
					listBCs, listClusterNames);
			
			if(iFile==0)
				totalp=p;
			else
				totalp.addPileup(p);
		}
		return totalp;
	}

	
	/**
	 * Build a pileup, everything in one group
	 */
	public Pileup buildPileup(
			String windowSeq, int windowFrom, int windowTo, 
			int numdiv) throws IOException {
		
		ArrayList<String> arrAllBC=new ArrayList<String>();
		ArrayList<String> arrAllFile=new ArrayList<String>();
		for(int i=0;i<listFile.size();i++) {
			CellPileFile cp=listFile.get(i);
			String pname=listPileName.get(i);
			
			List<String> list=cp.getListBarcodes();
			arrAllBC.addAll(list);
			for(int j=0;j<list.size();j++)
				arrAllFile.add(pname);
		}
		
		String[] cellBC=arrAllBC.toArray(new String[0]);
		String[] cellFile=arrAllFile.toArray(new String[0]);
		String[] cellCluster=getRepString("", cellBC.length);
		
		return buildPileup(windowSeq, windowFrom, windowTo, numdiv, cellBC, cellFile, cellCluster);
	}
	
	
	/**
	 * Return a string, repeated
	 */
	private static String[] getRepString(String el, int n) {
		String[] list=new String[n];
		for(int i=0;i<n;i++)
			list[i] = el;
		return list;
	}
	
}
