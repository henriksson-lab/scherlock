package isoform.cellpile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import isoform.counter.SplitGtf;

public class CellpileMain {
	
	public static void main(String[] args) throws IOException {
		
		
		
		File fCellpile=new File("/home/mahogny/temp/cellpile");
		File fBAM=new File("/beagle/big/henriksson/tonsil/possorted_genome_bam.bam");
		File fChromSizes=new File("/home/mahogny/umeå/project/isoform/refgenome/hg38.chrom.sizes");
		File fBC=new File("/beagle/big/henriksson/tonsil/filtered_feature_bc_matrix/barcodes.tsv.gz");
		
		ArrayList<String> listBarcodes=SplitGtf.readBarcodeZipList(fBC);
		
		CellpileFile cp=CellpileFile.writeFile(fCellpile, fChromSizes, fBAM, listBarcodes);
		cp.close();
		
		System.out.println("File written! ------ ");
		
		
		cp=CellpileFile.open(fCellpile);
		
		int[][] clusters=cp.convertBarcodeNamesToIDs(
				new String[][] {
						cp.getListBarcodes().toArray(new String[0])
				});
		
		int[][] tracks=cp.buildPileup(
				"1", 1, 100000, 
				100, 
				clusters);
		
		for(int[] track:tracks) {
			for(int i:track)
				System.out.print(i+" ");
			
			System.out.println();
		}
		
		
	}

}
