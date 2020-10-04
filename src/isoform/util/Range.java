package isoform.util;
/**
 * A range object, representing an exon or similar
 * 
 * @author Johan Henriksson
 *
 */
public class Range {
	public String source;
	public String featureType;
	public int from, to; //inclusive
	
	public String getSource() {return source;}
	public int getFrom() {return from;}
	public int getTo() {return to;}
	
	public Range(String source, int from, int to, String featureType) {
		this.source = source;
		this.from = from;
		this.to = to;
		this.featureType=featureType;
	}
	
	
	public final static String FEATURE_EXON="exon";
	public final static String FEATURE_5UTR="five_prime_UTR";
	public final static String FEATURE_3UTR="three_prime_UTR";

	/**
	 * Check if the feature overlaps a given range
	 */
	public boolean overlapsPos(int from, int to) {
		if(this.to<from || to<this.from)
			return false;
		else
			return true;
	}

	public boolean overlaps(int from, int to, String seq) {
		return this.source.equals(seq) && overlapsPos(from, to);
	}

	
	public boolean isUTR() {
		return featureType.equals(FEATURE_3UTR) || featureType.equals(FEATURE_5UTR);
	}
}
