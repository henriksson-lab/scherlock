package isoform.counter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.TreeMap;

import isoform.util.GtfParser;
import isoform.util.Range;


/**
 *
 * Split GTF file into suitable features
 * 
 * @author Johan Henriksson
 *
 */
public class GtfToFeature {
	
	/**
	 * Split GTF file into countable features
	 */
	public static ArrayList<Feature> splitGTF(File fGTF) throws IOException {
		ArrayList<Feature> features=new ArrayList<Feature>();
		TreeMap<String, ArrayList<Range>> mapGenesFeatures=new TreeMap<String, ArrayList<Range>>();

		GtfParser gtf=new GtfParser(fGTF);
		
		//Collapse transcript features into corresponding genes
		System.out.println("Merge transcripts");
		for(String tname:gtf.mapTranscriptGene.keySet()) {
			//Allocate gene array
			String gene=gtf.mapTranscriptGene.get(tname);
			ArrayList<Range> ga=mapGenesFeatures.get(gene);
			if(ga==null) {
				ga=new ArrayList<Range>();
				mapGenesFeatures.put(gene, ga);
			}
			
			//Transfer features from the transcript
			ga.addAll(gtf.mapTranscriptRanges.get(tname));
		}
		
		
		//Split features suitably, to handle casette exons etc
		System.out.println(gtf.mapGeneRange.size());
		System.out.println("Chopping up regions");
		for(String gname:new ArrayList<String>(mapGenesFeatures.keySet())) {
			
			//Sort the regions of a gene. Ordered first by from-pos, then to-pos
			ArrayList<Range> ta=mapGenesFeatures.get(gname);
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
						Range newr=new Range(r.source, lastTo+1, r.to, "");
						newta.add(newr);
						lastTo=r.to;					
					}
				}

				//Save it - list of region format ---  exon-like
				String source=ta.get(0).source;
				for(int i=0;i<newta.size();i++) {
					Range r=newta.get(i);
					String fname=gname+"_"+i+"_e";
					Feature cr=new Feature(gname, fname, "exon", source, r.from, r.to);
					features.add(cr);
				}
				
				//Save it - list of region format  --- intron-like
				for(int i=0;i<newta.size()-1;i++) {
					Range ra=newta.get(i);
					Range rb=newta.get(i);
					if(ra.to+1!=ra.from) {
						String fname=gname+"_"+i+"_i";
						Feature cr=new Feature(gname, fname, "intron",source, ra.to+1, rb.from-1);
						features.add(cr);
					}
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

		return features;
	}
	
	
	
	
	

}

