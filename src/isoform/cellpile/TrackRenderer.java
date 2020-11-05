package isoform.cellpile;

import java.util.ArrayList;

/**
 * 
 * One pileup, ready to be rendered
 * 
 * @author Johan Henriksson
 *
 */
public class TrackRenderer {
	public String seq;
	public int from;
	public int to; 
	
	public int numdiv;  //How many subdivisions the pileup was calculated at, from from-to
	
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
	double transformPosX(int x) {
		return labelsWidth + (double)trackWidth*(x-from)/(to-from);
	}
	
	double totalWidth=labelsWidth+trackWidth;

	ArrayList<Track> trackRenderer=new ArrayList<Track>();

	
	public boolean doLog;
	
	/**
	 * Render all of the pileup
	 */
	public String toSVG(boolean doLog) {
		this.doLog=doLog;
		
		StringBuilder sb=new StringBuilder();
		
		//Figure out total size
		double totalHeightAll=0;
		for(Track t:trackRenderer)
			totalHeightAll+=t.getHeight(this);//numTranscripts*transcriptHeight;
		
		//Write SVG header, with size
		sb.append("<svg height=\""+totalHeightAll+"\" width=\""+(totalWidth)+"\" "
				+ "style=\"background-color:none\">");  //seems to ignore bg
		
		//Mask for the track viewport
		sb.append(
				"<defs><mask id=\"trackmask\">" + 
				"<rect x=\""+labelsWidth+"\" y=\"0\" width=\""+totalWidth+"\" height=\""+totalHeightAll+"\" fill=\"white\"/>" + 
				"</mask></defs>");
		
		//Set background color
		sb.append("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>");
		
		renderVerticalGuideLines(sb);
		
		//Render all tracks
		int offsetY=0;
		for(Track t:trackRenderer) {
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
	 * Add the counts from another pileup. Assumes exactly the same settings etc or a crash will occur
	 */
	public void addPileup(TrackRenderer p) {
		for(int i=0;i<tracks.length;i++) {
			int[] t1=tracks[i];
			int[] t2=p.tracks[i];
			for(int j=0;j<t1.length;j++) {
				t1[j]+=t2[j];
			}
			clusterCellCount[i] += p.clusterCellCount[i];
		}
	}



	/**
	 * Add one track to render
	 */
	public void addTrack(TrackGTF t) {
		trackRenderer.add(t);
	}

}
