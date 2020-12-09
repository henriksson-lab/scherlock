package isoform.counter;

import java.util.TreeMap;

/**
 * A feature to be counted (genomic genome range)
 * 
 * @author Johan Henriksson
 *
 */
public class Feature {
	public String featureName;
	public String source;
	public int from, to; //inclusive
	
	TreeMap<String, String> extraAttr=new TreeMap<String, String>();

	
	
	private static String pop(TreeMap<String, String> a, String key) {
		if(a.containsKey(key)) {
			String val=a.get(key);
			a.remove(key);
			return val;
		} else {
			throw new RuntimeException("Missing feature key: "+key);
		}
	}
	
	

	public Feature(TreeMap<String, String> a) {
		featureName=pop(extraAttr,"feature");
		source=pop(extraAttr,"source");
		from=Integer.parseInt(pop(extraAttr,"from"));
		to=Integer.parseInt(pop(extraAttr,"to"));
		extraAttr.putAll(a);
	}

	
	public Feature(String gene, String fname, String type, String source, int from, int to) {
		this.featureName=fname;
		this.source = source;
		this.from = from;
		this.to = to;
		extraAttr.put("gene",gene);
		extraAttr.put("type",type);
	}

}
