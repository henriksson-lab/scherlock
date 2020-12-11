package isoform.counter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import isoform.util.PileUtil;


/**
 * Isoform counting, main class
 * 
 * @author Johan Henriksson
 *
 */
public class IsoCounterMain {

	/**
	 * Main function
	 */
	public static void main(String[] args) throws IOException {
		
		
		if(args.length ==3 && args[0].equals("build")) {		
//		if(args.length > 0 && args[0].equals("build")) {
//			if(args.length==3) {
				//////////////////////////////////////////////////////////////////
				////////////////// Build a basic list of features that can be counted
				//////////////////////////////////////////////////////////////////
				
				File fGTF=new File(args[1]);
				File fOut=new File(args[2]);
				
				if(!fGTF.exists()) {
					System.out.println("GTF-file does not exist");
					System.exit(1);
				}
			
				System.out.println("Building features from: "+fGTF);
				FeatureFile ff=FeatureFile.splitGTF(fGTF);
				
				System.out.println("Storing in: "+fOut);
				ff.write(fOut);
				
//			} else {
//				System.out.println("Too few arguments");
//			}
		} else if(args.length == 5 && args[0].equals("count")) {
//		} else if(args.length == 4 && args[0].equals("count")) {
//		} else if(args.length > 0 && args[0].equals("count")) {
//			if(args.length==4) {
				//////////////////////////////////////////////////////////////////
				////////////////// Perform counting
				//////////////////////////////////////////////////////////////////
				
				File fOut=new File(args[1]);
				File fFeature=new File(args[2]);
				File fBAM=new File(args[3]);
				File fBarcodes=new File(args[3]);
				
				if(!fFeature.exists()) {
					System.out.println("Feature-file does not exist: " + fFeature);
					System.exit(1);
				}

//				if(!f10x.exists()) {
//					System.out.println("10x directory does not exist: "+f10x);
//					System.exit(1);
//				}

				if(!fBAM.exists()) {
					System.out.println("10x directory does not exist: " + fBAM);
					System.exit(1);
				}
				
				if(!fBarcodes.exists()) {
					System.out.println("10x directory does not exist: " + fBarcodes);
					System.exit(1);
				}
				
				if(!fOut.exists()) {
					System.out.println("Outdir does not exist, creating: " + fOut);
					fOut.mkdir();
				}

//				//The regular 10x paths to filtered lists. Should we support unfiltered too?
//				File fBarcodes=new File(f10x, "filtered_feature_bc_matrix/barcodes.tsv.gz");
//				File fBAM=new File(f10x, "possorted_genome_bam.bam");

				System.out.println("Reading barcodes: "+fBarcodes);
				ArrayList<String> listBarcodes=PileUtil.readBarcodeZipList(fBarcodes);

				System.out.println("Reading features: "+fFeature);
				FeatureFile ff=FeatureFile.read(fFeature);
				
				System.out.println("Performing the counting: "+fBAM);
				CountFeaturesBAM cb=new CountFeaturesBAM(ff.features, listBarcodes);
				cb.countReads(fBAM);
				
				System.out.println("Storing the counts: "+fOut);
				cb.writeMatrix(fOut);			
				ff.writeZip(new File(fOut,"features.ext.tsv.gz"));
				
				
//			} else {
//				System.out.println("Too few arguments");
//			}
		} else { 
			System.out.println("==================================================");
			System.out.println("USAGE: isocounter build GTFFILE FEATUREFILE.ff");
			System.out.println("Where ");
			System.out.println("   GTFFILE is a something like Homo_sapiens.GRCh38.101.chr.gff3");
			System.out.println("   FEATUREFILE is where the features will be stored");
			System.out.println();
			System.out.println();
			System.out.println("==================================================");
			System.out.println("USAGE: isocounter count OUTDIR FEATUREFILE.ff sorted.bam barcodes.tsv.gz");
			System.out.println("Where ");
			System.out.println("   OUTDIR is where compressed new matrices will be stored");
			System.out.println("   FEATUREFILE.ff is a list of features to count");
			System.out.println("   sorted.bam is a position sorted bam file");
			System.out.println("   barcodes.tsv.gz is a gzipped tsv file with the cell barcodes in the first column");
			System.out.println();
			System.exit(1);
		}
	
	
	
	}
	
}

