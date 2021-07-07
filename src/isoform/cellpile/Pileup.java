package isoform.cellpile;

/**
 * 
 * One pileup, ready to be rendered
 * 
 * @author Johan Henriksson and Anton Björk
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
	public int[] clusterCellCount;  // This becomes slightly ambiguous since
	public String[] clusterNames;   // different cells can be represented among
									// the inbetweens and alignment blocks.
	//Pileups for each cluster 		// Left as is for now, ie not aware of inbetweens //AB
	public int[][] alignmentBlockTracks;
	public int[][] inbetweenTracks;
	
	
	
	/**
	 * Add the counts from another pileup. Assumes exactly the same settings etc or a crash will occur
	 */
	public void addPileup(Pileup p) {
		for(int i=0;i<alignmentBlockTracks.length;i++) {
			int[] t1=alignmentBlockTracks[i];
			int[] t2=p.alignmentBlockTracks[i];
			for(int j=0;j<t1.length;j++) {
				t1[j]+=t2[j];
			}
			clusterCellCount[i] += p.clusterCellCount[i];
		}

		for(int i=0;i<inbetweenTracks.length;i++) {
			int[] t1=inbetweenTracks[i];
			int[] t2=p.inbetweenTracks[i];
			for(int j=0;j<t1.length;j++) {
				t1[j]+=t2[j];
			}
		}
	}
}
