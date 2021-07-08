package isoform.trackRenderer;

import isoform.cellpile.Pileup;
import isoform.util.PileUtil;

/**
 * Renderer: Pileup track
 * 
 * @author Johan Henriksson and Anton Björk
 *
 */
public class TrackPileup extends Track {

	public Pileup pileup;
	public boolean showLog;
	
	
	
	public TrackPileup(Pileup pileup, boolean showLog) {
		this.pileup = pileup;
		this.showLog = showLog;
		this.trackName = "pileup";
	}

	@Override
	protected double getHeight(TrackRenderer renderer) {
		int numPileTrack=pileup.alignmentBlockTracks.length;
		return (int)(renderer.trackHeight*numPileTrack);
	}

	@Override
	protected void render(TrackRenderer renderer, StringBuilder sb, int offsetY) {
		
		//Rescale tracks
		double[][] tracksScaled=new double[pileup.alignmentBlockTracks.length][];
		double[] maxHeight=new double[pileup.alignmentBlockTracks.length];
		for(int curTrack=0;curTrack<pileup.alignmentBlockTracks.length;curTrack++) {
			//Rescale tracks vs #cells in track
			int[] thisTrack=pileup.alignmentBlockTracks[curTrack];
			double[] thisTrackScaled=tracksScaled[curTrack]=new double[thisTrack.length];
			double cnt=Math.max(1,pileup.clusterCellCount[curTrack]);
			for(int i=0;i<thisTrack.length;i++) {
				thisTrackScaled[i] = thisTrack[i]*cnt;
			}
			
			//Apply pseudo-log if requested
			if(showLog) {
				for(int i=0;i<thisTrack.length;i++) {
					thisTrackScaled[i] = Math.log10(1+thisTrackScaled[i]);
				}
			}
			maxHeight[curTrack]=PileUtil.maxForList(tracksScaled[curTrack]);
		}		
		double maxMaxHeight=Math.max(0.0000001, PileUtil.maxForList(maxHeight));

		//Write all pileup tracks
		for(int curTrack=0;curTrack<pileup.alignmentBlockTracks.length;curTrack++) {
			double baseX=renderer.labelsWidth;
			double baseY=(double)renderer.trackHeight*(curTrack+1);
			//scaleFactor[curTrack]/Math.max(1,cellGroupCellCount[curTrack]);
			///cellGroupCellCount[curTrack];//*maxForTrack[curTrack]);
			//*cellGroupCellCount[curTrack];
			// ?? Here before changes. //AB
			double scaleY=-(double)renderer.trackHeight*0.9/maxMaxHeight;
			double scaleX=(double)renderer.trackWidth/pileup.numdiv;
			renderOnePileup(renderer, sb, tracksScaled[curTrack], pileup.clusterNames[curTrack], baseX, baseY, scaleY, scaleX);
		}
	}

	
	

	
	/**
	 * Render one pileup track
	 */
	private void renderOnePileup(
			TrackRenderer renderer,
			StringBuilder sb,
			double[] track,
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
		double lastCount=-1;
		for(int i=0;i<track.length;i++) {
			double curCount=track[i];
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
		sb.append("<text x=\""+textXFrom+"\" y=\""+textY+"\" style=\""+textStyle+"\"  font-size=\""+renderer.textHeight+"px\" >"+trackName+"</text>");
	}

	protected void allocateSize(TrackRenderer renderer) {
	}
	
	
	
	
	
}
