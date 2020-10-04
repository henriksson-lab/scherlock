package isoform.util;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

/**
 *
 * Split GTF file into suitable features
 * 
 * @author Johan Henriksson
 *
 */
public class GtfParser {

	public TreeMap<String, Range> mapGeneRange=new TreeMap<String, Range>();
	public TreeMap<String, ArrayList<Range>> mapTranscriptRanges=new TreeMap<String, ArrayList<Range>>();

	public TreeMap<String,String> mapTranscriptGene=new TreeMap<String, String>();
	public TreeMap<String,Set<String>> mapGeneTranscripts=new TreeMap<String, Set<String>>();

	public TreeMap<String,String> mapSymbolGene=new TreeMap<String, String>();
	
	
	int numRange=0;
	
	/**
	 * Get list of transcripts for a gene.
	 * raw array for easy R/python interfacing
	 */
	public String[] getTranscriptsForGene(String gene) {
		return mapGeneTranscripts.get(gene).toArray(new String[0]);
	}
	
	/**
	 * Create a link gene <-> transcript.
	 * raw array for easy R/python interfacing
	 */
	private void linkGeneTranscript(String gene, String transcript) {
		mapTranscriptGene.put(transcript, gene);
		Set<String> m=mapGeneTranscripts.get(gene);
		if(m==null) {
			mapGeneTranscripts.put(gene,m=new TreeSet<String>());
		}
		m.add(transcript);
	}
	
	
	/**
	 * Get an array of [from to] for all the features of a transcript.
	 * raw array for easy R/python interfacing
	 */
	public int[] getTranscriptRangesFromTo(String transcript) {
		ArrayList<Range> list=mapTranscriptRanges.get(transcript);
		int[] outlist=new int[list.size()*2];
		for(int i=0;i<list.size();i++) {
			outlist[i*2+0]=list.get(i).from;
			outlist[i*2+1]=list.get(i).to;
		}
		return outlist;
	}
	
	/**
	 * Get an array of [featureType] for all the features of a transcript.
	 * raw array for easy R/python interfacing
	 */
	public String[] getTranscriptRangesFeatureType(String transcript) {
		ArrayList<Range> list=mapTranscriptRanges.get(transcript);
		String[] outlist=new String[list.size()];
		for(int i=0;i<list.size();i++) {
			outlist[i]=list.get(i).featureType;
		}
		return outlist;
	}


	/**
	 * Get the range of a gene. Gene can be symbol or ID
	 */
	public Range getRangeForGene(String gene) {
		if(mapSymbolGene.containsKey(gene))
			gene=mapSymbolGene.get(gene);
		return mapGeneRange.get(gene);
	}

	
	
	/**
	 * Get/Create array for transcripts
	 */
	private ArrayList<Range> getTranscriptArray(String gene){
		ArrayList<Range> a=mapTranscriptRanges.get(gene);
		if(a==null) {
			a=new ArrayList<Range>();
			mapTranscriptRanges.put(gene, a);
		}
		return a;
	}
	

	/**
	 * Parser with no entries
	 */
	public GtfParser() {
		
	}

	/**
	 * Parse out overall locations of all the genes
	 */
	public GtfParser(File fGTF) throws IOException {
		System.out.println("Reading features from "+fGTF);

		Reader reader;
		if(fGTF.getName().endsWith(".gz")) {
			reader=new InputStreamReader(new GZIPInputStream(new FileInputStream(fGTF)));
		} else {
			reader=new FileReader(fGTF);
		}
		
		//Read all the relevant features and associate with transcripts/genes
		BufferedReader br=new BufferedReader(reader);
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

				Range r=new Range(seq, sFrom, sTo, featureType);

				//System.out.println(attr);
				
				//Top level features
				String id=attr.get("ID");
				String parent=attr.get("Parent");
				if(id!=null) {
					if(id.startsWith("gene:")) {
						String attrGene=removeTag(id,"gene:");
						mapGeneRange.put(attrGene,r);
						if(!mapGeneTranscripts.containsKey(attrGene))
							mapGeneTranscripts.put(attrGene, new TreeSet<String>());
						
						String attrSymbol=attr.get("Name");
						if(attrSymbol!=null) {
							mapSymbolGene.put(attrSymbol, attrGene);
						}
						
					}

					if(id.startsWith("transcript:") && parent!=null && parent.startsWith("gene:")) {
						//TODO, filter
						String attrGene=removeTag(parent,"gene:");
						String attrTranscript=removeTag(id,"transcript:");
						linkGeneTranscript(attrGene, attrTranscript);
/*
						ArrayList<Range> ga=getTranscriptArray(attrTranscript);
						ga.add(r);
						numRange++;*/
					}
				}
				
				if(parent!=null) {
					//Exons, UTRs, etc
					if(parent.startsWith("transcript:")) {
						
						//Exons, UTRs, etc
						
						
						//TODO, filter
						String attrTranscript=removeTag(parent,"transcript:");
						
						ArrayList<Range> ga=getTranscriptArray(attrTranscript);
						ga.add(r);
						numRange++;
					}
					
				}
				
				
				
			}
		}
		br.close();		
		
		System.out.println("Stored ranges "+numRange);
		System.out.println("#genes "+mapGeneTranscripts.size());
		System.out.println("#transcript "+mapTranscriptGene.size());
		//Some genes might have no transcripts? check later, possibly special code in renderer
		
		//fixUTRoverlap();
	}
	
	/*
	private void fixUTRoverlap() {
		for(ArrayList<Range> trans:mapTranscriptRanges.values()) {
			if(!trans.isEmpty()) {
				
				if(trans.get(0).featureType.equals(Range.FEATURE_5UTR)) {
					
					
				}
				
				
			}
		}
	}*/

	
	
	
	/**
	 * Remove tag: in tag:FOO
	 */
	private String removeTag(String s, String tag) {
		if(s.startsWith(tag)) {
			return s.substring(tag.length());
		} else {
			throw new RuntimeException("Wrong tag; "+tag +" vs " +s);
		}
		
	}

	/**
	 * Parse an attribute string into a map
	 */
	private TreeMap<String,String> parseAttr(String s){
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

}
