package isoform.trackRenderer;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.util.Interval;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.bed.BEDFeature;

/**
 * Renderer: BED file track
 * 
 * @author Anton Björk, based on Johan Henrikssons TrackGTF.java
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
		
		// Had to add try catch to avoid syntax error (?) from unmatching/uncaught errors
		// Is there nicer way to do this?  //AB
		// Johan: No
		Interval interval = new Interval(renderer.seq, renderer.from, renderer.to);
		
		try {
			Iterator<BEDFeature> bedFileIterator;
			bedFileIterator = TrackBed.bedParser(BgzippedBedFileName, interval);
			// Make list of bed features/trackLines
			bedLineList = new ArrayList<TrackLine>();
			while (bedFileIterator.hasNext()) {
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

	
	
	public void render(TrackRenderer pileup, StringBuilder sb, int offsetY) {
//		String colorFeatureLine="rgb(0,0,0)";

		//Write BED annotation
		String textStyle="font: italic sans-serif; fill: red;";
		
		for (int ii = 0; ii<nTracks; ii++) {

//			ArrayList<BEDFeature> lineList = bedLineList.get(ii).featureList;
			TrackLine line = bedLineList.get(ii);
			String annotation = line.annotation;	
			BEDFeature feature = line.featureList.get(0);  // Could be a loop in future to go through several features sharing a track.

			double rFrom = pileup.transformPosX(feature.getStart());
			double rTo   = pileup.transformPosX(feature.getEnd());
			
//			//// Plot all the features
//			for(Range r:listFeatures) {
				//Transform coordinates
//				double rFrom=pileup.transformPosX(r.from);
//				double rTo=pileup.transformPosX(r.to);

			double lineY=offsetY + (ii+0.5)*pileup.transcriptHeight;
			double featureY=lineY - pileup.featureHeight/2;
			
			String exonColor="rgb(0,0,0)";
//				if(r.featureType.equals(Range.FEATURE_3UTR))
//					exonColor="rgb(0,255,0)";
//				if(r.featureType.equals(Range.FEATURE_5UTR))
//					exonColor="rgb(0,0,255)";
				
			String exonStyle="\"fill:"+exonColor+";stroke-width:3;stroke:none\"";

			sb.append("<rect x=\""+rFrom+"\" y=\""+featureY+"\" "
					+ "width=\""+(rTo-rFrom)+"\" "
					+ "height=\""+pileup.featureHeight+"\" style="+exonStyle+"/>");   //   mask=\"url(#trackmask)\"
				
//			}
			
			//Add text
			double textXFrom=5;
			double textY=lineY + pileup.textHeight/2;
			sb.append("<text x=\""+textXFrom+"\" y=\""+textY+"\" style=\""+textStyle+"\"  font-size=\""+pileup.textHeight+"px\" >"+annotation+"</text>");
		}
		
//			for(String gene:overlappingGenes) {
//				for(String transcript:gtf.mapGeneTranscripts.get(gene)) {
//					//System.out.println(transcript);
//					//Should draw a straight line...
//
//					double lineY=offsetY + (curTrack+0.5)*pileup.transcriptHeight;
//					double featureY=lineY - pileup.featureHeight/2;
//							
//					//// Plot the line also suggesting the direction
//					int minFrom=Integer.MAX_VALUE;
//					int maxTo=0;
//					for(Range r:gtf.mapTranscriptRanges.get(transcript)) {
//						minFrom=Math.min(minFrom,r.from);
//						maxTo=Math.max(maxTo,r.to);
//					}
////					double minTransFrom=transformPos(minFrom);
//		//			double maxTransTo=transformPos(maxTo);
//					double minTransFrom=Math.max(pileup.labelsWidth,pileup.transformPosX(minFrom));
//					double maxTransTo=Math.min(pileup.getWidth(),pileup.transformPosX(maxTo));
//
//					sb.append("<line"
//							+ " x1=\""+minTransFrom+
//							"\" y1=\""+lineY+
//							"\" x2=\""+maxTransTo+
//							"\" y2=\""+lineY+
//							"\" stroke=\""+
//							colorFeatureLine+"\" stroke-width=\"2px\"/>");  //mask=\"url(#trackmask)\"
//
//					//Plot the fishbones
//					
//					//// Prepare plotting: sort UTRs last as they overlap exons
//					ArrayList<Range> listFeatures=new ArrayList<Range>(gtf.mapTranscriptRanges.get(transcript));
//					listFeatures.sort(new Comparator<Range>() {
//						public int compare(Range a, Range b) {
//							boolean aUTR=a.isUTR();
//							boolean bUTR=b.isUTR();
//							if(aUTR) {
//								if(bUTR) {
//									return Integer.compare(a.from, b.from);
//								} else {
//									return 1;
//								}
//							} else {
//								if(bUTR) {
//									return -1;
//								} else {
//									return Integer.compare(a.from, b.from);
//								}
//							}
//						}
//					});
//					
//					//// Plot all the features
//					for(Range r:listFeatures) {
//						//Transform coordinates
//						double rFrom=pileup.transformPosX(r.from);
//						double rTo=pileup.transformPosX(r.to);
//
//						String exonColor="rgb(0,0,0)";
//						if(r.featureType.equals(Range.FEATURE_3UTR))
//							exonColor="rgb(0,255,0)";
//						if(r.featureType.equals(Range.FEATURE_5UTR))
//							exonColor="rgb(0,0,255)";
//						
//						String exonStyle="\"fill:"+exonColor+";stroke-width:3;stroke:none\"";
//						
//						sb.append("<rect x=\""+rFrom+"\" y=\""+featureY+"\" "
//								+ "width=\""+(rTo-rFrom)+"\" "
//								+ "height=\""+pileup.featureHeight+"\" style="+exonStyle+"/>");   //   mask=\"url(#trackmask)\"
//						
//						
//						//TO-DO style should depend on exon, utr
//					}
//
//					//TO-DO: add a line, and direction
//					
//				//Add text: Transcript ID
//				double textXFrom=5;
//				double textY=lineY + pileup.textHeight/2;
//				sb.append("<text x=\""+textXFrom+"\" y=\""+textY+"\" style=\""+textStyle+"\"  font-size=\""+pileup.textHeight+"px\" >"+transcript+"</text>");
//
//				curTrack++;
//				}
//			}
//	
			
			
		
		
	}

	@Override
	protected double getHeight(TrackRenderer renderer) {
		return nTracks*renderer.transcriptHeight;
	}

	
	
}
