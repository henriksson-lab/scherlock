package isoform.counter;

/**
 * A feature to be counted
 * 
 * @author Johan Henriksson
 *
 */
public class Feature {
	public String gene;
	public String featureName;
	
	public String source;
	public int from, to; //inclusive
	
	public Feature(String gene, String fname, String source, int from, int to) {
		this.gene=gene;
		this.featureName=fname;
		this.source = source;
		this.from = from;
		this.to = to;
	}

}
