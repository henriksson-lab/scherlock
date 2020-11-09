package isoform.trackRenderer;

import java.util.ArrayList;

/**
 * 
 * Renderer of genomic tracks
 * 
 * @author Johan Henriksson
 *
 */
public class TrackRenderer {
	
	//Where the camera should be
	public String seq;
	public int from;
	public int to; 
	
	//Rendering settings
	public double trackHeight=100;
	public double trackWidth=700;
	public double transcriptHeight=30;
	public double textHeight=20;
	public double featureHeight=15;
	public double labelsWidth=200;
	
	public int numVertGuides=16;
	
	public ArrayList<Track> tracks=new ArrayList<Track>();

	
	/**
	 * Create a renderer that immediately sets the view to correspond to the pileup
	 */
	public TrackRenderer(TrackPileup track) {
		addTrack(track);
		this.seq=track.pileup.seq;
		this.from=track.pileup.from;
		this.to=track.pileup.to;
	}

	/**
	 * Create a renderer, without a pileup track
	 */
	public TrackRenderer(String seq, int from, int to) {
		this.seq=seq;
		this.from=from;
		this.to=to;
	}
	
	/**
	 * Add one track to render
	 */
	public void addTrack(Track t) {
		tracks.add(t);
	}


	/**
	 * Transform position on chromosome to screen/SVG coordinate
	 */
	double transformPosX(int x) {
		return labelsWidth + (double)trackWidth*(x-from)/(to-from);
	}
	
	/**
	 * Total width of SVG
	 */
	public double getWidth() {
		return labelsWidth+trackWidth;
	}

	/**
	 * Total height of SVG
	 */
	public double getHeight() {
		double totalHeightAll=0;
		for(Track t:tracks)
			totalHeightAll+=t.getHeight(this);//numTranscripts*transcriptHeight;
		return totalHeightAll;
	}
	

	
	//public boolean doLog;
	
	/**
	 * Render all of the pileup
	 */
	public String toSVG() {
		
		//Allocate the size of all tracks
		for(Track t:tracks) {
			t.allocateSize(this);
		}
		
		//Write SVG header, with size
		double height=getHeight();
		double width=getWidth();
		StringBuilder sb=new StringBuilder();
		
		// Old line on master git branch
//		sb.append("<svg height=\""+height+"\" width=\""+width+"\" "
//				+ "style=\"background-color:none\">");  //seems to ignore bg
		// Merged
		sb.append("<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" height=\""+height+"\" width=\""+width+"\" "
				+ "style=\"background-color:none\">");  //seems to ignore bg
		// From fixed_svg_header git branch
//        sb.append("<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" height=\""
//                + totalHeightAll+"\" width=\""+(totalTrackWidth)+"\" "
//                + "style=\"background-color:none\">");  //seems to ignore bg

		
		
		//Mask for the track viewport
		sb.append(
				"<defs><mask id=\"trackmask\">" + 
				"<rect x=\""+labelsWidth+"\" y=\"0\" width=\""+width+"\" height=\""+height+"\" fill=\"white\"/>" + 
				"</mask></defs>");
		
		//Set background color
		sb.append("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>");
		
		renderVerticalGuideLines(sb);
		
		//Render all tracks
		int offsetY=0;
		for(Track t:tracks) {
			t.render(this, sb, offsetY);
			offsetY+=t.getHeight(this);
		}
		
		//End of file
		sb.append("</svg>");
		return sb.toString();
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
	 * Set all pileups to display log on/off
	 */
	public void setShowLog(boolean b) {
		for(Track t:tracks) {
			if(t instanceof TrackPileup) {
				((TrackPileup)t).showLog=b;
			}
		}
	}

	public void setTextHeight(double newTextHeight) {
		this.textHeight = newTextHeight;
	}

	
}

