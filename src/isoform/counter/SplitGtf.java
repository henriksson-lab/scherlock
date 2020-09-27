package isoform.counter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;


/**
 *
 * Split GTF file into suitable features
 * 
 * @author Johan Henriksson
 *
 */
public class SplitGtf {

	public TreeMap<String,String> mapTranscriptGene=new TreeMap<String, String>();
	public TreeMap<String, ArrayList<Range>> mapGenes=new TreeMap<String, ArrayList<Range>>();
	public TreeMap<String, ArrayList<Range>> mapTranscript=new TreeMap<String, ArrayList<Range>>();

	public ArrayList<Feature> features=new ArrayList<Feature>();
	
	
	public ArrayList<Range> getGeneArray(String gene){
		ArrayList<Range> a=mapGenes.get(gene);
		if(a==null) {
			a=new ArrayList<Range>();
			mapGenes.put(gene, a);
		}
		return a;
	}

	
	public ArrayList<Range> getTranscriptArray(String gene){
		ArrayList<Range> a=mapTranscript.get(gene);
		if(a==null) {
			a=new ArrayList<Range>();
			mapTranscript.put(gene, a);
		}
		return a;
	}

	
	/**
	 * Split GTF file into countable features
	 */
	public ArrayList<Feature> splitGTF(File fGTF, File fOutdir) throws IOException {
		System.out.println("Reading features from "+fGTF);

		//Read all the relevant features and associate with transcripts/genes
		BufferedReader br=new BufferedReader(new FileReader(fGTF));
		String line;
		while((line=br.readLine())!=null) {
			if(!line.startsWith("#")) {
				StringTokenizer stok=new StringTokenizer(line,"\t");
				
				String seq=stok.nextToken();
				stok.nextToken();
				String featureType=stok.nextToken();
				int sFrom=Integer.parseInt(stok.nextToken());
				int sTo=Integer.parseInt(stok.nextToken());
				stok.nextToken();
				stok.nextToken();//strand
				stok.nextToken();
				
				TreeMap<String,String> attr=parseAttr(stok.nextToken());

				Range r=new Range(seq, sFrom, sTo);

				if(featureType.contentEquals("gene")) {
					//Nothing to be done
				} else if(featureType.contentEquals("mRNA") || featureType.contentEquals("lnc_RNA")) {
					String attrTranscript=removeTag(attr.get("ID"),"transcript:");
					String attrGene=removeTag(attr.get("Parent"),"gene:");
					mapTranscriptGene.put(attrTranscript, attrGene);
				} else if(featureType.contentEquals("exon") || 
						featureType.contentEquals("three_prime_UTR") || 
						featureType.contentEquals("five_prime_UTR")) {
					//System.out.println(featureType);
					String attrParent=removeTag(attr.get("Parent"),"transcript:");
					ArrayList<Range> ga=getTranscriptArray(attrParent);
					ga.add(r);
				} /*else if(featureType.contentEquals("lnc_RNA")) {
					//System.out.println(featureType);
					String attrParent=removeTag(attr.get("Parent"),"gene:");
					ArrayList<Range> ga=getGeneArray(attrParent);
					ga.add(r);
				}*/else if(featureType.contentEquals("ncRNA_gene")) {
					String attrTranscript=removeTag(attr.get("ID"),"gene:");
					ArrayList<Range> ga=getGeneArray(attrTranscript);
					ga.add(r);
				}
				
			}
		}
		br.close();

		//Put transcript features into corresponding genes
		System.out.println("# transcript "+mapTranscriptGene.size());
		System.out.println("Merge transcripts");
		for(String tname:mapTranscriptGene.keySet()) {
			ArrayList<Range> ga=getGeneArray(mapTranscriptGene.get(tname));
			ArrayList<Range> ta=getTranscriptArray(tname);
			ga.addAll(ta);
		}
		
		//Split features suitably, to handle casette exons etc
		System.out.println(mapGenes.size());
		System.out.println("Chop up regions");
		for(String gname:new ArrayList<String>(mapGenes.keySet())) {
			
			//Sort the regions of a gene. Ordered first by from-pos, then to-pos
			ArrayList<Range> ta=getGeneArray(gname);
			ta.sort(new Comparator<Range>() {
				public int compare(Range ra, Range rb) {
					
					if(ra.from<rb.from) {
						return -1;
					} else if(ra.from>rb.from) {
						return 1;
					} else {
						return Integer.compare(ra.to, rb.to);
					}
				}
			});
			
			//Now do the splitting
			if(ta.isEmpty()) {
				System.out.println("Missing gene "+gname);
			} else {
				ArrayList<Range> newta=new ArrayList<Range>();
				int lastTo=ta.get(0).from-1;
				for(Range r:ta) {
					if(r.from>lastTo) {
						//Can add this range unmodified
						newta.add(r);
						lastTo=r.to;
					} else if (r.to <= lastTo) {
						//This part is covered; can skip entirely
					} else {
						//Need to add a chopped up version
						Range newr=new Range(r.source, lastTo+1, r.to);
						newta.add(newr);
						lastTo=r.to;					
					}
				}

				//Save it - list of region format
				String source=ta.get(0).source;
				int n=0;
				for(Range r:newta) {
					String fname=gname+"_"+n;
					Feature cr=new Feature(gname, fname, source, r.from, r.to);
					features.add(cr);
					n++;
				}
			}
		}		

		//Sort the features. This helps counting later on
		System.out.println("Position-sort features");
		features.sort(new Comparator<Feature>() {
			public int compare(Feature a, Feature b) {
				int sc = a.source.compareTo(b.source);
				if(sc==0) {
					return Integer.compare(a.from, b.from);
				} else {
					return sc;
				}
			}
		});

		//Store the features for visualization
		PrintWriter pwFeatures=new PrintWriter(new File(fOutdir,"features.bed"));
		for(Feature r:features) {
			pwFeatures.println(r.source+"\t"+r.from+"\t"+r.to+"\t"+r.featureName);
			
		}
		pwFeatures.close();

		
		return features;
		
	}
	
	
	/**
	 * Remove tag: in tag:FOO
	 */
	public String removeTag(String s, String tag) {
		if(s.startsWith(tag)) {
			return s.substring(tag.length());
		} else {
			throw new RuntimeException("Wrong tag; "+tag +" vs " +s);
		}
		
	}

	/**
	 * Parse an attribute string into a map
	 */
	public TreeMap<String,String> parseAttr(String s){
		TreeMap<String,String> m=new TreeMap<String, String>();
		StringTokenizer stok=new StringTokenizer(s,";");
		while(stok.hasMoreTokens()) {
			String tok=stok.nextToken();
			
			int ind=tok.indexOf("=");
			String key=tok.substring(0,ind);
			String value=tok.substring(ind+1);
			m.put(key, value);
		}
		return m;
	}
	
	
	/**
	 * Read list of zip file e.g. barcode list
	 */
	public static ArrayList<String> readBarcodeZipList(File fBarcodes) throws IOException{
		ArrayList<String> list=new ArrayList<String>();
		BufferedReader brFeatures=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fBarcodes))));
		String line;
		while((line=brFeatures.readLine())!=null) {
			list.add(line);
		}
		brFeatures.close();
		return list;
	}
	
	
	
	
	
	
	

}

