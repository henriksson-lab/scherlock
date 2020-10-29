package isoform.cellpile;


import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.bed.BEDCodec.StartOffset;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class TestReadBedFile {
	public static void main(String[] args) throws Exception {  // Why do I have to declare exceptions??
		
		// Read first byte of file from input stream
//		int smallRead = is.read();
//		System.out.println(smallRead);
		
		// Read and print contents of bed file to console
//		AsciiLineReader lr = new AsciiLineReader(is);  // Java docs: Use readline() method instead. I think they mean the one in BufferedReader
		String filePath = "C:\\Users\\anton\\java_projects\\isocounter\\data\\UP000005640_9606_lipid.bed";
		
		// I can use a FileReader instead of these two
//		InputStream is = new java.io.FileInputStream(filePath);
//		InputStreamReader isr = new InputStreamReader(is);
		
		// Parse BED file
		List<BEDFeature> allFeatures;
		allFeatures = new ArrayList<BEDFeature>();
		HashMap<String, List<BEDFeature>> chrToBF;
		chrToBF = new HashMap<String, List<BEDFeature>>();
	    BEDCodec bc = new BEDCodec(StartOffset.ZERO);
		FileReader fr = new FileReader(filePath);
		BufferedReader br = new BufferedReader(fr);
		String cl = null;
		while((cl = br.readLine()) != null) {
			BEDFeature bf = bc.decode(cl);
			if(bf != null)
			{
				// List of all features
				allFeatures.add(bf);
				// HashMap from contig to feature.
//            	if(chrToBF.containsKey(bf.getChr()))  // getChr deprecated for getContig
				if(chrToBF.containsKey(bf.getContig()))  
				{
					chrToBF.get(bf.getContig()).add(bf); 
				}
				else
				{
					List<BEDFeature> a = new ArrayList<BEDFeature>();
					a.add(bf);
					chrToBF.put(bf.getContig(), a);
				}
			}
		}
		
		// Print chr9 keyset
		System.out.println(chrToBF.keySet());
		
		// Print names of all features in chr9
		List<BEDFeature> chrom = chrToBF.get("chr9");
		for (int ii = 0; ii<chrom.size(); ii++) {
			System.out.println(chrom.get(ii).getLink());
		}
		
		fr.close();
		br.close();
	}
}



