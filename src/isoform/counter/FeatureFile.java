package isoform.counter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.StringTokenizer;


/**
 * 
 * File to keep track of regions in which to count reads
 * 
 * @author Johan Henriksson
 *
 */
public class FeatureFile {
	
	
	public ArrayList<Feature> features=new ArrayList<Feature>();

	/**
	 * Generate a file from a GTF
	 */
	public static FeatureFile splitGTF(File fGTF) throws IOException {
		FeatureFile ff=new FeatureFile();
		ff.features=GtfToFeature.splitGTF(fGTF);
		return ff;
	}
	
	/**
	 * Read a feature file
	 */
	public static FeatureFile read(File fFeature) throws IOException {
		FeatureFile ff=new FeatureFile();
		BufferedReader br=new BufferedReader(new FileReader(fFeature));
		String line;
		while((line=br.readLine())!=null) {
			StringTokenizer stok=new StringTokenizer(line, "\t");
			String gene=stok.nextToken();
			String fname=stok.nextToken();
			String source=stok.nextToken();
			int from=Integer.parseInt(stok.nextToken());
			int to=Integer.parseInt(stok.nextToken());
			
			Feature f=new Feature(gene, fname, source, from, to);
			ff.features.add(f);
		}
		br.close();
		return ff;
	}


	/**
	 * Store the features
	 */
	public void write(File outfile) throws IOException {
		PrintWriter pwFeatures=new PrintWriter(outfile);
		for(Feature r:features) {
			pwFeatures.println(r.gene+"\t"+r.featureName+"\t"+r.source+"\t"+r.from+"\t"+r.to);
		}
		pwFeatures.close();
	}
	
	

}
