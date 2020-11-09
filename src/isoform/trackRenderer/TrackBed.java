package isoform.trackRenderer;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import EDU.oswego.cs.dl.util.concurrent.FJTask.Seq;
import htsjdk.samtools.util.Interval;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.bed.BEDFeature;

/**
 * Renderer: BED file track
 * 
 * @author Anton Bj�rk, based on Johan Henrikssons TrackGTF.java
 *
 */
public class TrackBed extends Track {

	public TrackBed(String BgzippedBedFileName) {
		this.BgzippedBedFileName = BgzippedBedFileName;
	}

	
	String BgzippedBedFileName;

	
	
	//Allocation: which genes to show?
	//	private int numTranscripts;  					rename --> nTracks 
	//	private TreeSet<String> overlappingGenes;  		rename --> bedLineList 
	private int nTracks;
	private List<TrackLine> bedLineList;
	
	
	
	private static class TrackLine {
		ArrayList<BEDFeature> featureList = new ArrayList<BEDFeature>();
		String annotation;
	}
	
	
//	public static void  main(String[] args) throws Exception {
//		
////		String bedBgzipFileName = "merged_all.sorted.bed_details_compatible.bed.gz";
////		
////		// cd55 location
////		Interval interval = new Interval("chr1", 207321376, 207360966);
//		
//		
//	}


	//
//	public TrackBed(GtfParser gtf) {
//		this.gtf = gtf;
//		this.trackName = "gtf";
//	}
	
	// For bed file parsing //AB
//	public GtfParser gtf;
	public static Iterator<BEDFeature> bedParser(String BgzippedBedFileName, Interval interval) throws Exception{
		final FeatureReader<BEDFeature> reader = AbstractFeatureReader.getFeatureReader(BgzippedBedFileName, new BEDDetailCodec());
		final Iterator<BEDFeature> readerIterator = 
				reader.query(interval.getContig(), interval.getStart(), interval.getEnd());
		System.out.println(interval.getContig());
		System.out.println(interval.getStart());
		System.out.println(interval.getEnd());
		
//		while (readerIterator.hasNext()) {
//			BEDFeature bedFeature = readerIterator.next();
////			System.out.println(bedFeature);
////			System.out.println(bedFeature.getColor());
////			System.out.println(bedFeature.getDescription());
////			System.out.println(bedFeature.getExons());
////			System.out.println(bedFeature.getLink());
////			System.out.println(bedFeature.getName());
////			System.out.println(bedFeature.getScore());
////			System.out.println(bedFeature.getStrand());
////			System.out.println(bedFeature.getType());
////			System.out.println(bedFeature.getStart());
////			System.out.println(bedFeature.getEnd());
////			System.out.println();	
//		}
		
		return readerIterator;
	}

		
	
	/**
	 * Get BED features to show
	 * Return void but sets nTracks and bedTrackList
	 */
	public void allocateSize(TrackRenderer renderer) {
		
		System.out.println("Entered allocateSize method");
		// Had to add try catch to avoid syntax error (?) from unmatching/uncaught errors
		// Is there nicer way to do this?  //AB
		// Johan: No
		
		// Temp translation bc incompatible test data. To be removed after fixing test data. //AB
//		final String chrom = "chr" + renderer.seq;
//		System.out.println("chrom");
//		Interval interval = new Interval(chrom, renderer.from, renderer.to);
		
		Interval interval = new Interval(renderer.seq, renderer.from, renderer.to);
		
		try {  // This try/catch block is to prevent Java from complaining about uncaught errors
			 Iterator<BEDFeature> bedFileIterator = TrackBed.bedParser(BgzippedBedFileName, interval);
			// Make list of bed features/trackLines
			bedLineList = new ArrayList<TrackLine>();
			while (bedFileIterator.hasNext()) {
				System.out.println("Hurra!");
				BEDFeature bedFeature = bedFileIterator.next();

				TrackLine trackLine	= new TrackLine();
				trackLine.annotation = bedFeature.getName();  // We might want to change this to .getDescription() later,
															  // 	unless we can get description as a hover pop-up			
				// If we want several features on the same track later, 
				// we can do something nice here to add relevant ones together on the same line. 
				// Hence a list of a single feature for now. Type wont change if we implement that change later. "Future-proof"
				ArrayList<BEDFeature> bedFeatureList = new ArrayList<BEDFeature>();
				bedFeatureList.add(bedFeature);
				trackLine.featureList = bedFeatureList;
				
				bedLineList.add(trackLine);
			}	
			nTracks = bedLineList.size();
			System.out.println(nTracks);
			System.out.println(bedLineList);
		} 
		catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}

		
	}

	
	
	public void render(TrackRenderer renderer, StringBuilder sb, int offsetY) {

		String textStyle="font: italic sans-serif; fill: red;";
		int canvasWidthBp = renderer.to - renderer.from;
		
		
		for (int ii = 0; ii<nTracks; ii++) {

			TrackLine line = bedLineList.get(ii);
			String annotation = line.annotation;	
			BEDFeature feature = line.featureList.get(0);  // Could be a loop in future to go through several features sharing a track.

			// Rescale feature if too small to see.
			int rFromRaw = feature.getStart();
			int rToRaw   = feature.getEnd();
			int minimum_feature_size = canvasWidthBp/1000;
			if (rToRaw - rFromRaw < minimum_feature_size) {
				rToRaw = rFromRaw + minimum_feature_size;  // Otherwise, feature will be almost invisible in plot
			}
			double rFrom = renderer.transformPosX(rFromRaw);
			double rTo   = renderer.transformPosX(rToRaw);

			// Set y positions
			double lineY=offsetY + (ii+0.5)*renderer.transcriptHeight;
			double featureY=lineY - renderer.featureHeight/2;
			
			// Add rectangle
			String rectColor="rgb(0,0,0)";
			String rectStyle="\"fill:"+rectColor+";stroke-width:3;stroke:none\"";
			sb.append("<rect x=\""+rFrom+"\" y=\""+featureY+"\" "
					+ "width=\""+(rTo-rFrom)+"\" "
					+ "height=\""+renderer.featureHeight+"\" style="+rectStyle+"/>");   //   mask=\"url(#trackmask)\"
			
			//Add text
			double textXFrom=5;
			double textY=lineY + renderer.textHeight/2;
			sb.append("<text x=\""+textXFrom+"\" y=\""+textY+"\" style=\""+textStyle+"\"  font-size=\""+renderer.textHeight+"px\" >"+annotation+"</text>");
		}
	}

	@Override
	protected double getHeight(TrackRenderer renderer) {
		return nTracks*renderer.transcriptHeight;
	}

	
	
}
