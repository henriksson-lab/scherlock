package isoform.cellpile;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import htsjdk.samtools.util.Interval;
import isoform.trackRenderer.TrackBed;
import isoform.trackRenderer.TrackGTF;
import isoform.trackRenderer.TrackPileup;
import isoform.trackRenderer.TrackRenderer;
import isoform.util.GtfParser;

/**
 * Test reading functions
 * 
 * @author Johan Henriksson
 *
 */
public class CellPileTestRead {
	
	public static void main(String[] args) throws IOException {
		
//		// Parse cellpile
//		File fCellpile=new File("/home/mahogny/temp/cellpile");
//		CellPileFile cp=CellPileFile.open(fCellpile);
//		int[][] clusters=cp.convertBarcodeNamesToIDs(
//				new String[][] {
//						cp.getListBarcodes().toArray(new String[0])
//				});
//		String[] clusterNames=new String[] {"all"};
//		//cd55: "1",207321376,207360966
//		Pileup pileup=cp.buildPileup(
//				"1",207321376,207360966,
//				//"1", 1, 50000, 
//				1000, 
//				clusters, clusterNames);

		
//		long time=System.currentTimeMillis();
		
//		// Parse GTF
//		//File fGTF=new File("/home/mahogny/umeå/project/isoform/refgenome/Homo_sapiens.GRCh38.101.chr.gff3");
//		File fGTF=new File("/home/mahogny/umeå/project/isoform/refgenome/Homo_sapiens.GRCh38.101.chr.gff3.gz");
////		File fGTF=new File("/home/mahogny/all.gtf.gz");
//		System.out.println("parse gtf");
//		GtfParser gtf=new GtfParser(fGTF);
////		System.out.println("dt: "+(System.currentTimeMillis()-time));
//		
//		// Render pileup
//		System.out.println("--------");
//		TrackPileup trackPileup=new TrackPileup(pileup, false);
//		TrackRenderer renderer=new TrackRenderer(trackPileup);
//		
//		// Render GTF
//		renderer.addTrack(new TrackGTF(gtf));
		
		
		String BgzippedBedFileName = "C:\\Users\\anton\\java_projects\\isocounter\\src\\isoform\\trackRenderer\\merged_all.sorted.bed_details_compatible.bed.gz";
//		 cd55 location
		Interval interval = new Interval("chr1", 207321376, 207360966);
		String chromosome = interval.getContig();
		int start = interval.getStart();
		int end = interval.getEnd();
		
		
		// Make empty render
		TrackRenderer renderer=new TrackRenderer(chromosome, start, end);
		
		// Render BED file
		renderer.addTrack(new TrackBed(BgzippedBedFileName));
		
		
		// Convert to svg
		String svg=renderer.toSVG();
		//System.out.println(svg);
		
		// Print svg
		System.out.println("--------");
		PrintWriter pw=new PrintWriter(new File("test.svg"));
		pw.print(svg);
		pw.close();
		
		
		// Done
		System.out.println("done");
		
		/*
		for(int[] track:pileup.tracks) {
			for(int i:track)
				System.out.print(i+" ");
			System.out.println();
		}*/
		
		
	}

}
