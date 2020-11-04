package isoform.cellpile;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.TreeSet;

import isoform.util.GtfParser;
import isoform.util.Range;

/**
 * 
 * One pileup, ready to be rendered
 * 
 * @author Johan Henriksson
 *
 */
public class Pileup {
	public String seq;
	public int from;
	public int to; 
	public int numdiv;
	public int[][] cellGroups;
	public int[] cellGroupCellCount;
	public int[][] tracks;
	public String[] clusterNames;
	
	
	//Rendering settings
	public double trackHeight=100;
	public double trackWidth=700;
	public double transcriptHeight=30;
	public double textHeight=20;
	public double featureHeight=15;
	public double labelsWidth=200;
	
	public int numVertGuides=16;
	
	/**
	 * Transform position on chromosome to screen/SVG coordinate
	 */
	private double transformPosX(int x) {
		return labelsWidth + (double)trackWidth*(x-from)/(to-from);
	}
	
	
	/**
	 * Render one pileup track
	 */
	private void renderOneTrackToSVG(
			StringBuilder sb,
			int[] track,
			String trackName,
			double baseX,
			double baseY,
			double scaleY,
			double scaleX
			) {
		
		String stylePolygon="fill:black; stroke:black;";
		sb.append("<polygon points=\"");
		
		//Move to the first point
		sb.append(""+baseX+","+baseY+" ");
		
		boolean needToAddLast=false;
		int lastCount=-1;
		for(int i=0;i<track.length;i++) {
			int curCount=track[i];
			if(curCount!=lastCount) {
				if(needToAddLast)
					sb.append(""+(baseX+(i-1)*scaleX)+","+(baseY+lastCount*scaleY)+" ");
				
				//Only add a new point if different from the last. Speeds up rendering
				sb.append(""+(baseX+i*scaleX)+","+(baseY+curCount*scaleY)+" ");
				lastCount=curCount;
				needToAddLast=false;
			} else {
				needToAddLast=true;
			}
		}
		//All the way to the right side
		sb.append(""+(baseX+track.length*scaleX)+" "+(baseY+track[track.length-1]*scaleY)+" ");
		//To lower right corner
		sb.append(""+(baseX+track.length*scaleX)+" "+(baseY)+" ");
		//Close and end
		sb.append("\" style=\""+stylePolygon+"\"/>"); 

		
		String textStyle="font: italic sans-serif; fill: red;";
		double textXFrom=5;
		double textY=baseY;
		sb.append("<text x=\""+textXFrom+"\" y=\""+textY+"\" style=\""+textStyle+"\"  font-size=\""+textHeight+"px\" >"+trackName+"</text>");
	}
	
	
	/**
	 * Render all of the pileup
	 */
	public String toSVG(GtfParser gtf) {
		StringBuilder sb=new StringBuilder();
		
		//Figure out the space needed for tracks
		int numTrack=tracks.length;

		//Figure out which genes and transcripts to show - linear search, not optimal
		TreeSet<String> overlappingGenes=new TreeSet<String>();
		int numTranscripts=0;
		for(String gene:gtf.mapGeneRange.keySet()) {
			Range range=gtf.mapGeneRange.get(gene);
			if(range.overlaps(from,to,seq)) {
				overlappingGenes.add(gene);
				numTranscripts+=gtf.mapGeneTranscripts.get(gene).size();
			}
		}
		//System.out.println("#transcript "+numTranscripts);

		//Write SVG header, with size
		double totalTrackWidth=labelsWidth+trackWidth;
		int totalTrackHeight=(int)(trackHeight*numTrack);
		double totalHeightAll=(totalTrackHeight + numTranscripts*transcriptHeight);
		sb.append("<svg height=\""+totalHeightAll+"\" width=\""+(totalTrackWidth)+"\" "
				+ "style=\"background-color:none\">");  //seems to ignore bg
		
		//Mask for the track viewport
		sb.append(
				"<defs><mask id=\"trackmask\">" + 
				"<rect x=\""+labelsWidth+"\" y=\"0\" width=\""+totalTrackWidth+"\" height=\""+totalHeightAll+"\" fill=\"white\"/>" + 
				"</mask></defs>");
		
		//Set background color
		sb.append("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>");
		
		renderVerticalGuideLines(sb);
		renderAllTracksToSVG(sb);
		renderGTF(sb, gtf, overlappingGenes, totalTrackWidth, totalTrackHeight);
		
		
		sb.append("</svg>");
		return sb.toString();
	}
	
	
	
	private void renderGTF(
			StringBuilder sb, 
			GtfParser gtf, TreeSet<String> overlappingGenes, 
			double totalTrackWidth, double totalTrackHeight) {
		String colorFeatureLine="rgb(0,0,0)";

		//Write GTF annotation
		String textStyle="font: italic sans-serif; fill: red;";
		int curTrack=0;
		for(String gene:overlappingGenes) {
			for(String transcript:gtf.mapGeneTranscripts.get(gene)) {
				//System.out.println(transcript);
				//Should draw a straight line...

				double lineY=totalTrackHeight + (curTrack+0.5)*transcriptHeight;
				double featureY=lineY - featureHeight/2;
						
				//// Plot the line also suggesting the direction
				int minFrom=Integer.MAX_VALUE;
				int maxTo=0;
				for(Range r:gtf.mapTranscriptRanges.get(transcript)) {
					minFrom=Math.min(minFrom,r.from);
					maxTo=Math.max(maxTo,r.to);
				}
//				double minTransFrom=transformPos(minFrom);
	//			double maxTransTo=transformPos(maxTo);
				double minTransFrom=Math.max(labelsWidth,transformPosX(minFrom));
				double maxTransTo=Math.min(totalTrackWidth,transformPosX(maxTo));

				sb.append("<line"
						+ " x1=\""+minTransFrom+
						"\" y1=\""+lineY+
						"\" x2=\""+maxTransTo+
						"\" y2=\""+lineY+
						"\" stroke=\""+
						colorFeatureLine+"\" stroke-width=\"2px\"/>");  //mask=\"url(#trackmask)\"

				//Plot the fishbones
				
				//// Prepare plotting: sort UTRs last as they overlap exons
				ArrayList<Range> listFeatures=new ArrayList<Range>(gtf.mapTranscriptRanges.get(transcript));
				listFeatures.sort(new Comparator<Range>() {
					public int compare(Range a, Range b) {
						boolean aUTR=a.isUTR();
						boolean bUTR=b.isUTR();
						if(aUTR) {
							if(bUTR) {
								return Integer.compare(a.from, b.from);
							} else {
								return 1;
							}
						} else {
							if(bUTR) {
								return -1;
							} else {
								return Integer.compare(a.from, b.from);
							}
						}
					}
				});
				
				//// Plot all the features
				for(Range r:listFeatures) {
					//Transform coordinates
					double rFrom=transformPosX(r.from);
					double rTo=transformPosX(r.to);

					String exonColor="rgb(0,0,0)";
					if(r.featureType.equals(Range.FEATURE_3UTR))
						exonColor="rgb(0,255,0)";
					if(r.featureType.equals(Range.FEATURE_5UTR))
						exonColor="rgb(0,0,255)";
					
					String exonStyle="\"fill:"+exonColor+";stroke-width:3;stroke:none\"";
					
					sb.append("<rect x=\""+rFrom+"\" y=\""+featureY+"\" "
							+ "width=\""+(rTo-rFrom)+"\" "
							+ "height=\""+featureHeight+"\" style="+exonStyle+"/>");   //   mask=\"url(#trackmask)\"
					
					
					//TODO style should depend on exon, utr
				}

				//TODO: add a line, and direction
				
			//Add text: Transcript ID
			double textXFrom=5;
			double textY=lineY + textHeight/2;
			sb.append("<text x=\""+textXFrom+"\" y=\""+textY+"\" style=\""+textStyle+"\"  font-size=\""+textHeight+"px\" >"+transcript+"</text>");

			curTrack++;
			}
		}
	}

	
	/**
	 * Render the vertical guide lines
	 */
	private void renderVerticalGuideLines(StringBuilder sb) {
		String colorVertGuides="rgb(200,200,200)";
		//String colorVertGuides="lightgray";
		//Write vertical guide lines
		for(int i=0;i<numVertGuides;i++) {
			double x=labelsWidth + i*trackWidth/numVertGuides;
			sb.append("<line x1=\""+x+"\" y1=\"0\" x2=\""+x+"\" y2=\"100%\" stroke=\""+colorVertGuides+"\"/>");
		}
	}

	/**
	 * Render all tracks. The main challenge is to scale them properly
	 */
	private void renderAllTracksToSVG(StringBuilder sb) {
		
		//Figure out track rescaling vs number of cells in each vs max peak height
		double[] scaleFactor=new double[tracks.length];
		double[] scaleFactorWithHeight=new double[tracks.length];
		double[] maxHeight=new double[tracks.length];
		for(int curTrack=0;curTrack<tracks.length;curTrack++) {
			System.out.println("cell count: "+cellGroupCellCount[curTrack]);
			maxHeight[curTrack]=Math.max(1, PileMathUtil.maxForList(tracks[curTrack]));
			scaleFactor[curTrack] = 1.0/Math.max(1,cellGroupCellCount[curTrack]);
			scaleFactorWithHeight[curTrack] = Math.max(1,cellGroupCellCount[curTrack])/maxHeight[curTrack];
			//scaleFactor[curTrack] /= Math.max(1.0,maxHeight[curTrack]);//Math.max(1,cellGroupCellCount[curTrack]);
		}
		PileMathUtil.scaleMaxValue(scaleFactor, PileMathUtil.maxForList(scaleFactorWithHeight));   //maybe min?
		/*
		//Figure out track rescaling vs max peak height, taking #cell rescaling into account
		double[] scaleByCellsAndMaximum=new double[tracks.length];
		for(int curTrack=0;curTrack<tracks.length;curTrack++) {
			scaleByCellsAndMaximum[curTrack]=scaleByCells[curTrack]/Math.max(1.0,maxForList(tracks[curTrack]));
		}
		scaleMaxValueToOne(scaleByCellsAndMaximum);*/
		
		//Write all pileup tracks
		for(int curTrack=0;curTrack<tracks.length;curTrack++) {
			double baseX=labelsWidth;
			double baseY=(double)trackHeight*(curTrack+1);
			double scaleY=-(double)trackHeight*scaleFactor[curTrack]/Math.max(1,cellGroupCellCount[curTrack]);///cellGroupCellCount[curTrack];//*maxForTrack[curTrack]);//*cellGroupCellCount[curTrack];
			double scaleX=(double)trackWidth/numdiv;
			renderOneTrackToSVG(sb, tracks[curTrack], clusterNames[curTrack], baseX, baseY, scaleY, scaleX);
		}
	}

	



//	public double maxForTrack[];
	
	/**
	 * Find the maximum count on any track
	 * 
	 * After normalization by 1/cellGroupCellCount[iTrack]   -- has to be undone later
	 * 
	 */
	/*
	private double getMaxCountForTracks() {
		double max=0;
		maxForTrack=new double[tracks.length];
		for(int iTrack=0;iTrack<tracks.length;iTrack++) {
			int[] track=tracks[iTrack];
			
			//Calculate normalized maximum for this track
			maxForTrack[iTrack]=maxForList(track);
			///(double)(Math.max(1,cellGroupCellCount[iTrack]));    //should normalize on total counts!
			System.out.println("count:    "+cellGroupCellCount[iTrack]);
			
			max=Math.max(maxForTrack[iTrack], max);
		}
		return max;
	}*/
	
	
	/**
	 * Add the counts from another pileup. Assumes exactly the same settings etc or a crash will occur
	 */
	public void addPileup(Pileup p) {
		for(int i=0;i<tracks.length;i++) {
			int[] t1=tracks[i];
			int[] t2=p.tracks[i];
			for(int j=0;j<t1.length;j++) {
				t1[j]+=t2[j];
			}
			
			cellGroupCellCount[i] += p.cellGroupCellCount[i];
		}
	}

}
