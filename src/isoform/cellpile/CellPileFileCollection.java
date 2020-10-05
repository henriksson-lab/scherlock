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
	 * Get the names of piles
	 */
	public List<String> getPileNames(){
		return listPileName;
	}
	
	
	/**
	 * Convert barcodes to IDs for a given fileID
	 */
	private int[] convertBarcodeNamesToIDs(
			int forPile, 
			String[] cellBC, String[] cellPile, String[] cellCluster, 
			String forCluster){
		CellPileFile pile=listFile.get(forPile);
		String thisPileName=listPileName.get(forPile);
		
		
		IntArrayList list2=new IntArrayList(cellBC.length);
		for(int j=0;j<cellBC.length;j++) {
			if(cellPile[j].equals(thisPileName) && cellCluster[j].equals(forCluster)) {
				list2.add(pile.mapBarcodeIndex.get(cellBC[j]));
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
