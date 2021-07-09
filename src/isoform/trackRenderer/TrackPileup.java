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
	// trackname not in here even though set by constructor. Problem or purpose? //AB
	public boolean showInbetweens;
	private int nAlignBlockTracks;
	private int nTracks;

	
	public TrackPileup(Pileup pileup, boolean showLog, boolean showInbetweens) {
		this.pileup = pileup;
		this.showLog = showLog;
		this.trackName = "pileup";
		this.showInbetweens = showInbetweens;
		this.nAlignBlockTracks = pileup.alignmentBlockTracks.length;
		// Always one inbetween track for every alignment block track
		this.nTracks = showInbetweens ? this.nAlignBlockTracks*2 : this.nAlignBlockTracks;
	}

	@Override
	protected double getHeight(TrackRenderer renderer) {
		return (int)(renderer.trackHeight*this.nTracks);
	}



	@Override
	protected void render(TrackRenderer renderer, StringBuilder sb, int offsetY) {

		double[][] tracksScaled=new double[this.nTracks][];
		double[] maxHeight=new double[this.nTracks];

		int[][] tracks;
		int[] nCellsTracks;
		String[] trackNames;

		if (showInbetweens) {
			tracks = reorderTracks();
			nCellsTracks = reorderNumberOfCellsPerTrack();
			trackNames = expandTrackNames();
		} else {
			tracks = this.pileup.alignmentBlockTracks;
			nCellsTracks = this.pileup.clusterCellCount;
			trackNames = this.pileup.clusterNames;
		}

		//Rescale tracks
		// for(int curTrack=0;curTrack<pileup.alignmentBlockTracks.length;curTrack++) {  // Old version
		for(int curTrack=0;curTrack<nTracks;curTrack++) {
			//Rescale tracks vs #cells in track
			int[] thisTrack=tracks[curTrack];
			double[] thisTrackScaled=tracksScaled[curTrack]=new double[thisTrack.length];
			double cnt=Math.max(1, nCellsTracks[curTrack]);
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
		for(int curTrack=0;curTrack<tracks.length;curTrack++) {
			double baseX=renderer.labelsWidth;
			double baseY=(double)renderer.trackHeight*(curTrack+1);
			
			//scaleFactor[curTrack]/Math.max(1,cellGroupCellCount[curTrack]);
			///cellGroupCellCount[curTrack];//*maxForTrack[curTrack]);
			//*cellGroupCellCount[curTrack];
			// ?? Here before changes. //AB

			double scaleY;
			double scaleX;
			if (renderer.individualTrackScaling) {
				scaleY=-(double)renderer.trackHeight*0.9/maxHeight[curTrack];
				scaleX=(double)renderer.trackWidth/pileup.numdiv;
			} else {
				scaleY=-(double)renderer.trackHeight*0.9/maxMaxHeight;
				scaleX=(double)renderer.trackWidth/pileup.numdiv;
			}
			renderOnePileup(renderer, sb, tracksScaled[curTrack], 
				trackNames[curTrack], baseX, baseY, scaleY, scaleX);
		}
	}


	// Reorders tracks so that alignmentblock tracks are
	// directly above corresponding inbetween tracks.
	private int[][] reorderTracks() {
		int[][] tracksReordered=new int[this.nTracks][];
		for(int curTrack=0;curTrack<this.nAlignBlockTracks;curTrack++) {
			tracksReordered[curTrack*2] = this.pileup.alignmentBlockTracks[curTrack];
			tracksReordered[curTrack*2+1] = this.pileup.inbetweenTracks[curTrack];
		}
		return tracksReordered;
	}

	// Expands the cell counts for the alignment block tracks to also cover
	// the inbetween tracks. Should be used with reorderTracks() since
	// the ordering of nCellsTracks follow that of tracksReordered
	private int[] reorderNumberOfCellsPerTrack() {
		int[] nCellsTracks=new int[this.nTracks];
		for(int curTrack=0;curTrack<this.nAlignBlockTracks;curTrack++) {
			nCellsTracks[curTrack*2] = this.pileup.clusterCellCount[curTrack];
			nCellsTracks[curTrack*2+1] = this.pileup.clusterCellCount[curTrack];
		}
		return nCellsTracks;
	}

	// Expands the cell counts for the alignment block tracks to also cover
	// the inbetween tracks. Should be used with reorderTracks() since
	// the ordering of nCellsTracks follow that of tracksReordered
	private String[] expandTrackNames() {
		String[] out=new String[this.nTracks];
		for(int curTrack=0;curTrack<this.nAlignBlockTracks;curTrack++) {
			out[curTrack*2] = this.pileup.clusterNames[curTrack] + " Aligned Blocks";
			out[curTrack*2+1] = this.pileup.clusterNames[curTrack] + " Inbetweens";
		}
		return out;
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
