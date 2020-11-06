package isoform.cellpile;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

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
		
		File fCellpile=new File("C:\\Users\\anton\\java_projects\\isocounter\\data\\SRR11816791.cellpile");
		
		CellPileFile cp=CellPileFile.open(fCellpile);
		
		int[][] clusters=cp.convertBarcodeNamesToIDs(
				new String[][] {
						cp.getListBarcodes().toArray(new String[0])
				});
		String[] clusterNames=new String[] {"all"};
		

		//cd55: "1",207321376,207360966
		
		
		Pileup pileup=cp.buildPileup(
				"1",207321376,207360966,
				//"1", 1, 50000, 
				1000, 
				clusters, clusterNames);

		long time=System.currentTimeMillis();
		//File fGTF=new File("/home/mahogny/umeå/project/isoform/refgenome/Homo_sapiens.GRCh38.101.chr.gff3");
		File fGTF=new File("C:\\Users\\anton\\java_projects\\isocounter\\data\\Homo_sapiens.GRCh38.101.chr.gff3.gz");
//		File fGTF=new File("/home/mahogny/all.gtf.gz");
		
		System.out.println("parse gtf");
		GtfParser gtf=new GtfParser(fGTF);
		System.out.println("dt: "+(System.currentTimeMillis()-time));
		
		System.out.println("--------");
		TrackPileup trackPileup=new TrackPileup(pileup, false);
		
		TrackRenderer renderer=new TrackRenderer(trackPileup);
		renderer.addTrack(new TrackGTF(gtf));
		String svg=renderer.toSVG();
		//System.out.println(svg);
		System.out.println("--------");
		PrintWriter pw=new PrintWriter(new File("out\\test.svg"));
		pw.print(svg);
		pw.close();
		System.out.println("done");
		
		/*
		for(int[] track:pileup.tracks) {
			for(int i:track)
				System.out.print(i+" ");
			System.out.println();
		}*/
		
		
	}

}
