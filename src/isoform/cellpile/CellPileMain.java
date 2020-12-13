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
		if(cmd.equals("build") && args.length==5) {

			File fCellpile=new File(args[1]);
			File fChromSizes=new File(args[2]);
			File fBAM=new File(args[3]);
			File fBC=new File(args[4]);
			
			ArrayList<String> listBarcodes=PileUtil.readBarcodeZipList(fBC);
			
			System.out.println("Building cellpile from: "+fBAM);
			System.out.println("To: "+fCellpile);
			CellPileFile cp=CellPileFile.writeFile(fCellpile, fChromSizes, fBAM, listBarcodes);
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
		} else {
			
			System.out.println("Usages: ");
			System.out.println("=============================================================================");
			System.out.println("java -jar cellpile.jar build OUTFILE.cellpile sizes.chromosome sorted.bam barcodes.tsv.gz");
			System.out.println("   where");
			System.out.println("      sorted.bam is a position sorted bam file");
			System.out.println("      barcodes.tsv.gz is a gzipped tsv file with the cell barcodes in the first column");
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
