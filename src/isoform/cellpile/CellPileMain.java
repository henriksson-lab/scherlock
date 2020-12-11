package isoform.cellpile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import isoform.util.PileUtil;

public class CellPileMain {
	
	public static void main(String[] args) throws IOException {
		
		String cmd="";
		if(args.length>0) {
			cmd=args[0];
		}
		if(cmd.equals("build") && args.length==5) {

			File fCellpile=new File(args[1]);
			File fChromSizes=new File(args[2]);
//			File f10x=new File(args[3]);
			
//			File fBAM=new File(f10x,"possorted_genome_bam.bam");
//			File fBC=new File(f10x, "filtered_feature_bc_matrix/barcodes.tsv.gz");
			
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
			
		} else {
			
			System.out.println("Usages: ");
			System.out.println("=============================================================================");
			System.out.println("java -jar cellpile.jar build OUTFILE.cellpile sizes.chromosome sorted.bam barcodes.tsv.gz");
			System.out.println("   where");
			System.out.println("   sorted.bam is a position sorted bam file");
			System.out.println("   barcodes.tsv.gz is a gzipped tsv file with the cell barcodes in the first column");
			System.out.println("");
			System.out.println("java -jar cellpile.jar inspect OUTFILE.cellpile");
			System.out.println("=============================================================================");
			System.out.println("How to produce chromosome sizes:");
			System.out.println("  samtools faidx input.fa");
			System.out.println("  cut -f1,2 input.fa.fai > sizes.genome");
			
		}
	}
}
