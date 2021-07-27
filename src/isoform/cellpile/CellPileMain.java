package isoform.cellpile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import isoform.util.PileUtil;
import isoform.util.RandomSubset;

public class CellPileMain {
	
	public static void main(String[] args) throws IOException {
		
		String cmd="";
		if(args.length>0) {
			cmd=args[0];
		}

		// Argument parsing
		// args[5] can be either File or String 
		// depending on option --single_cell vs --bulk
		// No matter, it is used to populate listBarcodes
		if(cmd.equals("build") && args.length==6 && args[1].equals("--single_cell")) {

			String bamType = args[1];
			File fCellpile=new File(args[2]);
			File fChromSizes=new File(args[3]);
			File fBam=new File(args[4]);
			File fBC=new File(args[5]);	
			
			ArrayList<String> listBarcodes=PileUtil.readBarcodeZipList(fBC);
			
			System.out.println("Building cellpile from: "+fBam);
			System.out.println("To: "+fCellpile);
			CellPileFile cp=CellPileFile.writeFile(fCellpile, fChromSizes, fBam, bamType, listBarcodes);
			cp.close();
			
			System.out.println("Done writing cellpile");

		} else if(cmd.equals("build") && args.length==6 && args[1].equals("--bulk")) {

			String bamType = args[1];
			File fCellpile=new File(args[2]);
			File fChromSizes=new File(args[3]);
			File fBam=new File(args[4]);
			String bulk_id=new String(args[5]);
			
			// New bulk version; only one barcode
			ArrayList<String> listBarcodes = new ArrayList<String>();
			listBarcodes.add(bulk_id);
			
			System.out.println("Building cellpile from: "+fBam);
			System.out.println("To: "+fCellpile);
			CellPileFile cp=CellPileFile.writeFile(fCellpile, fChromSizes, fBam, bamType, listBarcodes);
			cp.close();
			
			System.out.println("Done writing cellpile");

		} else if(cmd.equals("inspect") && args.length==2) {
			
			File fCellpile=new File(args[1]);
			CellPileFile cp=CellPileFile.open(fCellpile);

			System.out.println("Number of barcodes included: "+cp.getListBarcodes().size());
			System.out.println("Number of sequences included: "+cp.getListSequences().size());

		} else if(cmd.equals("randomsubset") && args.length==6) {

			File fBarcodes=new File(args[1]);
			File fInBam=new File(args[2]);
			File fOutBam=new File(args[3]);
			double pKeepBC=Double.parseDouble(args[4]);
			double pKeepNonBC=Double.parseDouble(args[5]);

			System.out.println("Reading barcodes: "+fBarcodes);
			ArrayList<String> listBarcodes=PileUtil.readBarcodeZipList(fBarcodes);

			System.out.println("Subsetting");
			RandomSubset.subset(fInBam, fOutBam, new HashSet<String>(listBarcodes), pKeepBC, pKeepNonBC);

		} else if(cmd.equals("pileinfo") && args.length==2) {

			File fInPile=new File(args[1]);

                        CellPileFile.printInfo(fInPile);

		} else {
			
			System.out.println("Usages: ");
			System.out.println("=============================================================================");
			System.out.println("java -jar cellpile.jar build bamType outfile.cellpile sizes.chromosome sorted.bam barcodes");
			System.out.println("   where");
			System.out.println("      bamType is either '--single_cell' for single cell RNAseq data, or '--bulk' for bulk RNAseq data");
			System.out.println("      sorted.bam is a position sorted bam file");
			System.out.println("      bamType:");
			System.out.println("          If bamType is 'single_cell', barcodes is a gzipped tsv file with the cell barcode whitelist in the first column");
			System.out.println("          If bamType is 'bulk', barcodes is a string with a name/barcode for the entire bulk sample");
			System.out.println("");
			System.out.println("=============================================================================");
			System.out.println("java -jar cellpile.jar inspect OUTFILE.cellpile");
			System.out.println("   How to produce chromosome sizes:");
			System.out.println("     samtools faidx input.fa");
			System.out.println("     cut -f1,2 input.fa.fai > sizes.genome");
			System.out.println("");
			System.out.println("=============================================================================");
			System.out.println("java -jar cellpile.jar randomsubset BARCODES.tsv.gz IN.bam OUT.bam P_KEEP_BC P_KEEP_NONBC");
			System.out.println("   Produces random subsamplings. P-values given as e.g. 0.05");
			System.out.println();
		}
	}
}
