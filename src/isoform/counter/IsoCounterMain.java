package isoform.counter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;


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
		
		//File fFeatures=new File("/beagle/big/debojyoti/HUMAN_SARS_229E_OC43/lib1/outs/filtered_feature_bc_matrix/features.tsv.gz");
		//ArrayList<String> listFeatures=readZipList(fFeatures);
		
		//File fGTF=new File("/home/mahogny/umeå/project/isoform/refgenome/Homo_sapiens.GRCh38.101.chr.gff3");
		//File f10x=new File("/beagle/big/henriksson/tonsil");
		//File fOut=new File("/home/mahogny/umeå/project/isoform/newcount");

		///home/mahogny/umeå/project/isoform/refgenome/Homo_sapiens.GRCh38.101.chr.gff3 /beagle/big/henriksson/tonsil /home/mahogny/umeå/project/isoform/newcount

		
	//		File f10x=new File("/beagle/big/debojyoti/HUMAN_SARS_229E_OC43/lib1/outs");
		if(args.length==4 && args[0].equals("build")) {
			File fGTF=new File(args[1]);
			File fOut=new File(args[2]);
			
			if(!fGTF.exists()) {
				System.out.println("GTF-file does not exist");
				System.exit(1);
			}

			System.out.println("Reading "+fGTF);
			FeatureFile ff=FeatureFile.read(fGTF);
			ff.write(fOut);
			System.out.println("Stored in "+fOut);
			
		} else if(args.length==4 && args[0].equals("count")) {

			File fFeature=new File(args[1]);
			File f10x=new File(args[2]);
			File fOut=new File(args[3]);

			if(!fFeature.exists()) {
				System.out.println("Feature-file does not exist");
				System.exit(1);
			}

			if(!f10x.exists()) {
				System.out.println("10x directory does not exist");
				System.exit(1);
			}

			if(!fOut.exists()) {
				System.out.println("Outdir does not exist");
				System.exit(1);
			}

			//scp rackham.uppmax.uu.se:/home/mahogny/mystore/dataset/tonsil_sc/aligned/SRR11816791/outs/possorted_genome_bam.bam
			File fBarcodes=new File(f10x, "filtered_feature_bc_matrix/barcodes.tsv.gz");
			File fBAM=new File(f10x, "possorted_genome_bam.bam");

			System.out.println("Reading barcodes: "+fBarcodes);
			ArrayList<String> listBarcodes=GtfToFeature.readBarcodeZipList(fBarcodes);

			System.out.println("Reading features: "+fFeature);
			FeatureFile ff=FeatureFile.read(fFeature);
			
			System.out.println("Performing the counting: "+fBAM);
			CountFeaturesBAM cb=new CountFeaturesBAM(ff.features, listBarcodes);
			cb.countReads(fBAM);
			
			System.out.println("Storing the counts: "+fOut);
			cb.writeMatrix(fOut);
			
		} else {
			System.out.println("==================================================");
			System.out.println("USAGE: isocounter build GTFFILE FEATUREFILE.ff");
			System.out.println("Where ");
			System.out.println("   GTFFILE is a something like Homo_sapiens.GRCh38.101.chr.gff3");
			System.out.println("   FEATUREFILE is where the features will be stored");
			System.out.println();
			System.out.println();
			System.out.println("==================================================");
			System.out.println("USAGE: isocounter count FEATUREFILE.ff MATRIXDIR OUTDIR");
			System.out.println("Where ");
			System.out.println("   FEATUREFILE is a list of features to count");
			System.out.println("   MATRIXDIR is the 10x matrix out-directory");
			System.out.println("   OUTDIR is where compressed new matrices will be stored, and what the features are (.bed)");
			System.out.println();
			System.exit(1);
		}
	
	
	
	}
	
}
