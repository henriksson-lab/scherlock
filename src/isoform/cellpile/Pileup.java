package isoform.cellpile;

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
	public int[][] tracks;
	public String[] clusterNames;
	
	
	//Rendering settings
	public double oneTrackH=100;
	public double oneTrackW=800;
	public double oneTranscriptH=30;
	public double textHeight=20;
	public double featureHeight=20;
	public double tracksBeginX=200;
	
	public int numVertGuides=16;
	
	private double transformPos(int x) {
		return tracksBeginX + (double)oneTrackW*(x-from)/(to-from);
	}
	
	private void trackToSVG(
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
		
		int lastCount=-1;
		for(int i=0;i<track.length;i++) {
			int curCount=track[i];
			if(curCount!=lastCount) {
				//Only add a new point if different from the last. Speeds up rendering
				sb.append(""+(baseX+i*scaleX)+","+(baseY+curCount*scaleY)+" ");
				lastCount=curCount;
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
		double textY=baseX;
		//+ (curt*oneTranscriptH) + featureHeight - (textHeight-featureHeight)/2;
		sb.append("<text x=\""+textXFrom+"\" y=\""+textY+"\" style=\""+textStyle+"\"  font-size=\""+textHeight+"px\" >"+trackName+"</text>");

		
	}
	
	
	public String toSVG(GtfParser gtf) {
		StringBuilder sb=new StringBuilder();

		
		//Figure out the space needed for tracks
		int numTrack=tracks.length;
		int trackMaxCount=Math.max(1,getMaxCountForTracks());

		//Figure out which genes to show - linear search, not optimal
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
		double totalTrackWidth=tracksBeginX+oneTrackW;
		int totalTrackHeight=(int)(oneTrackH*numTrack);
		double totalHeightAll=(totalTrackHeight + numTranscripts*oneTranscriptH);
		sb.append("<svg height=\""+totalHeightAll+"\" width=\""+(totalTrackWidth)+"\" "
				+ "style=\"background-color:none\">");  //seems to ignore bg
		
		sb.append("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>");
		
		
		String colorVertGuides="rgb(200,200,200)";
		//String colorVertGuides="lightgray";
		//Write vertical guide lines
		for(int i=0;i<numVertGuides;i++) {
			double x=tracksBeginX + i*oneTrackW/numVertGuides;
			sb.append("<line x1=\""+x+"\" y1=\"0\" x2=\""+x+"\" y2=\"100%\" stroke=\""+colorVertGuides+"\"/>");
		}
		
		
		//Write all tracks
		for(int curTrack=0;curTrack<tracks.length;curTrack++) {
			double baseX=tracksBeginX;
			double baseY=(double)oneTrackH*(curTrack+1);
			double scaleY=-(double)oneTrackH/trackMaxCount;
			double scaleX=(double)oneTrackW/numdiv;
			trackToSVG(sb, tracks[curTrack], clusterNames[curTrack], baseX, baseY, scaleY, scaleX);
		}
		
		
		//Write GTF annotation
		String exonStyle="\"fill:rgb(0,0,0);stroke-width:3;stroke:none\"";
		String textStyle="font: italic sans-serif; fill: red;";
		int curt=0;
		for(String gene:overlappingGenes) {
			for(String transcript:gtf.mapGeneTranscripts.get(gene)) {
				//System.out.println(transcript);
				//Should draw a straight line...
				
				for(Range r:gtf.mapTranscriptRanges.get(transcript)) {
					//Transform coordinates
					double rFrom=transformPos(r.from);
					double rTo=transformPos(r.to);
					double y=totalTrackHeight + curt*oneTranscriptH + (oneTranscriptH-featureHeight)/2;

					
					
					sb.append("<rect x=\""+rFrom+"\" y=\""+y+"\" "
							+ "width=\""+(rTo-rFrom)+"\" "
							+ "height=\""+featureHeight+"\" style="+exonStyle+" />");
					
					
					//TODO style should depend on exon, utr
				}

				//TODO: add a line, and direction
				
			//Add text: Transcript ID
			double textXFrom=5;
			double textY=totalTrackHeight + (curt*oneTranscriptH) + featureHeight - (textHeight-featureHeight)/2;
			sb.append("<text x=\""+textXFrom+"\" y=\""+textY+"\" style=\""+textStyle+"\"  font-size=\""+textHeight+"px\" >"+transcript+"</text>");
				

			curt++;
			}
		}

		
		
		
		sb.append("</svg>");
		return sb.toString();
	}
	
	/**
	 * Find the maximum count on any track
	 */
	private int getMaxCountForTracks() {
		int max=0;
		for(int[] track:tracks) {
			for(int x:track) {
				max=Math.max(max, x);
			}
		}
		return max;
	}
	
	
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
		}
	}

}
