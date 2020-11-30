package isoform.counter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.StringTokenizer;
import java.util.zip.GZIPOutputStream;


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
			String source=stok.nextToken();
			int from=Integer.parseInt(stok.nextToken());
			int to=Integer.parseInt(stok.nextToken());
			String fname=stok.nextToken();
			String gene=stok.nextToken();
			String type=stok.nextToken();
			
			Feature f=new Feature(gene, fname, type, source, from, to);
			ff.features.add(f);
		}
		br.close();
		
		//Sort the entries
		ff.features.sort(new Comparator<Feature>() {
			public int compare(Feature a, Feature b) {
				int comp=a.source.compareTo(b.source);
				if(comp==0) {
					return Integer.compare(a.from, b.from);
				} else {
					return comp;
				}
			}
		});
		
		
		
		return ff;
	}


	/**
	 * Store the features
	 */
	public void write(File outfile) throws IOException {
		FileOutputStream os=new FileOutputStream(outfile);
		write(os);
		os.close();
	}
	

	/**
	 * Store the features
	 */
	public void write(OutputStream os) throws IOException {
		PrintWriter pwFeatures=new PrintWriter(os);
		for(Feature r:features) {
			pwFeatures.println(r.source+"\t"+r.from+"\t"+r.to+"\t"+r.featureName+"\t"+r.gene+"\t"+r.type);
		}
		pwFeatures.close();
	}

	public void writeZip(File file) throws IOException {
		OutputStream os=new GZIPOutputStream(new FileOutputStream(file));
		write(os);
		os.close();
	}
	

}
