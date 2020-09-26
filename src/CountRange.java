
public class CountRange {
	String gene;
	String featureName;
	
	String source;
	int from, to; //inclusive
	
	public CountRange(String gene, String fname, String source, int from, int to) {
		this.gene=gene;
		this.featureName=fname;
		this.source = source;
		this.from = from;
		this.to = to;
	}

}
