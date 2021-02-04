package isoform.cellpile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import isoform.util.PileUtil;

/**
 * Test functions
 * 
 * @author Johan Henriksson
 *
 */
public class CellPileTestWrite {
	
	public static void main(String[] args) throws IOException {
		
		
		String bamType="--single_cell";
		File fCellpile=new File("/home/mahogny/temp/cellpile");
		File fBAM=new File("/big/henriksson/tonsil/possorted_genome_bam.bam");
		File fChromSizes=new File("/data/henlab/ref_genome/human/sizes.genome");
		File fBC=new File("/big/henriksson/tonsil/filtered_feature_bc_matrix/barcodes.tsv.gz");
		
		ArrayList<String> listBarcodes=PileUtil.readBarcodeZipList(fBC);
		
		CellPileFile cp=CellPileFile.writeFile(fCellpile, fChromSizes, fBAM, bamType, listBarcodes);
		cp.close();
		
		System.out.println("File written! ------ ");
		
		/*
		cp=CellPileFile.open(fCellpile);
		
		int[][] clusters=cp.convertBarcodeNamesToIDs(
				new String[][] {
						cp.getListBarcodes().toArray(new String[0])
				});
		
		Pileup pileup=cp.buildPileup(
				"1", 1, 100000, 
				100, 
				clusters);
		
		System.out.println("");
		System.out.println(pileup.toSVG());
		*/
		/*
		for(int[] track:pileup.tracks) {
			for(int i:track)
				System.out.print(i+" ");
			
			System.out.println();
		}*/
		
		
	}

}
