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
import java.util.TreeMap;
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
		
		//Read header
		ArrayList<String> header=new ArrayList<String>();
		line=br.readLine();
		StringTokenizer stok2=new StringTokenizer(line, "\t");
		while(stok2.hasMoreTokens()) {
			header.add(stok2.nextToken());
		}
		
		//Read features
		while((line=br.readLine())!=null) {
			StringTokenizer stok=new StringTokenizer(line, "\t");			
			int i=0;
			TreeMap<String, String> map=new TreeMap<String, String>();
			while(stok.hasMoreTokens()) {
				map.put(header.get(i),stok.nextToken());
				i++;
			}

			Feature f=new Feature(map);
			ff.features.add(f);
		}
		br.close();
		
		ff.sortByPosition();
		return ff;
	}


	/**
	 * Sort entries by position, same order as sorted BAM files
	 */
	private void sortByPosition() {
		features.sort(new Comparator<Feature>() {
			public int compare(Feature a, Feature b) {
				int comp=a.source.compareTo(b.source);
				if(comp==0) {
					return Integer.compare(a.from, b.from);
				} else {
					return comp;
				}
			}
		});
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
		
		//Write header
		pwFeatures.print("source\tfrom\tto\tfeature");
		if(!features.isEmpty()) {
			for(String n:features.get(0).extraAttr.keySet()) {
				pwFeatures.print("\t"+n);
			}
		}
		pwFeatures.println();
		
		//Write all the features
		for(Feature r:features) {
			pwFeatures.println(r.source+"\t"+r.from+"\t"+r.to+"\t"+r.featureName);
			for(String val:r.extraAttr.values()) {
				pwFeatures.print("\t"+val);
			}
			pwFeatures.println();
		}
		pwFeatures.close();
	}

	public void writeZip(File file) throws IOException {
		OutputStream os=new GZIPOutputStream(new FileOutputStream(file));
		write(os);
		os.close();
	}
	

}
