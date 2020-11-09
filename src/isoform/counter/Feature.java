package isoform.counter;

/**
 * A feature to be counted (genomic genome range)
 * 
 * @author Johan Henriksson
 *
 */
public class Feature {
	public String gene;
	public String featureName;
	
	public String source;
	public String type;
	public int from, to; //inclusive
	
	public Feature(String gene, String fname, String type, String source, int from, int to) {
		this.gene=gene;
		this.featureName=fname;
		this.type=type;
		this.source = source;
		this.from = from;
		this.to = to;
	}

}
