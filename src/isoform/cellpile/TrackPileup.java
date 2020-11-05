package isoform.cellpile;

public class TrackPileup extends Track {

	@Override
	protected double getHeight(Pileup pileup) {
		int numPileTrack=pileup.tracks.length;
		return (int)(pileup.trackHeight*numPileTrack);
	}

	@Override
	protected void render(Pileup pileup, StringBuilder sb, int offsetY) {
		
		//Rescale tracks
		double[][] tracksScaled=new double[pileup.tracks.length][];
		double[] maxHeight=new double[pileup.tracks.length];
		for(int curTrack=0;curTrack<pileup.tracks.length;curTrack++) {
			//Rescale tracks vs #cells in track
			int[] thisTrack=pileup.tracks[curTrack];
			double[] thisTrackScaled=tracksScaled[curTrack]=new double[thisTrack.length];
			double cnt=Math.max(1,pileup.clusterCellCount[curTrack]);
			for(int i=0;i<thisTrack.length;i++) {
				thisTrackScaled[i] = thisTrack[i]*cnt;
			}
			
			//Apply pseudo-log if requested
			if(pileup.doLog) {
				for(int i=0;i<thisTrack.length;i++) {
					thisTrackScaled[i] = Math.log10(1+thisTrackScaled[i]);
				}
			}
			maxHeight[curTrack]=PileMathUtil.maxForList(tracksScaled[curTrack]);
		}		
		double maxMaxHeight=Math.max(0.0000001, PileMathUtil.maxForList(maxHeight));

		//Write all pileup tracks
		for(int curTrack=0;curTrack<pileup.tracks.length;curTrack++) {
			double baseX=pileup.labelsWidth;
			double baseY=(double)pileup.trackHeight*(curTrack+1);
			double scaleY=-(double)pileup.trackHeight*0.9/maxMaxHeight;//scaleFactor[curTrack]/Math.max(1,cellGroupCellCount[curTrack]);///cellGroupCellCount[curTrack];//*maxForTrack[curTrack]);//*cellGroupCellCount[curTrack];
			double scaleX=(double)pileup.trackWidth/pileup.numdiv;
			renderOnePileup(pileup, sb, tracksScaled[curTrack], pileup.clusterNames[curTrack], baseX, baseY, scaleY, scaleX);
		}
	}

	
	

	
	/**
	 * Render one pileup track
	 */
	private void renderOnePileup(
			Pileup pileup,
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
		sb.append("<text x=\""+textXFrom+"\" y=\""+textY+"\" style=\""+textStyle+"\"  font-size=\""+pileup.textHeight+"px\" >"+trackName+"</text>");
	}
	
	
	
	
	
}
