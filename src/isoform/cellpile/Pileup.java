package isoform.cellpile;

/**
 * 
 * One pileup, ready to be rendered
 * 
 * @author Johan Henriksson
 *
 */
public class Pileup {
	//Where the pileup is calculated
	public String seq;
	public int from;
	public int to; 
	
	public int numdiv;  //How many subdivisions the pileup was calculated at, from from-to
	
	//The cluster info
	public int[][] cellCluster;
	public int[] clusterCellCount;
	public String[] clusterNames;
	
	//Pileups for each cluster
	public int[][] tracks;
	
	
	
	/**
	 * Add the counts from another pileup. Assumes exactly the same settings etc or a crash will occur
	 */
	public void addPileup(Pileup p) {
		for(int i=0;i<tracks.length;i++) {
			int[] t1=tracks[i];
			int[] t2=p.tracks[i];
			for(int j=0;j<t1.length;j++) {
				t1[j]+=t2[j];
			}
			clusterCellCount[i] += p.clusterCellCount[i];
		}
	}



}
