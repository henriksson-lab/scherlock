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
	
	public Range(String source, int from, int to, String featureType) {
		this.source = source;
		this.from = from;
		this.to = to;
		this.featureType=featureType;
	}
	
	
	public final static String FEATURE_EXON="exon";
	public final static String FEATURE_5UTR="five_prime_UTR";
	public final static String FEATURE_3UTR="three_prime_UTR";

}
