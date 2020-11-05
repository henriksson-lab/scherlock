package isoform.cellpile;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.TreeSet;

import isoform.util.GtfParser;
import isoform.util.Range;

/**
 * Track that shows a GTF file
 * 
 * @author Johan Henriksson
 *
 */
public class TrackGTF extends Track {

	GtfParser gtf;
	


	public TrackGTF(GtfParser gtf) {
		this.gtf = gtf;
	}

	
	
	int numTranscripts;
	TreeSet<String> overlappingGenes;

	/**
	 * Figure out which genes and transcripts to show
	 */
	public void allocateSize(Pileup pileup) {
		numTranscripts=0;
				
		//linear search, not optimal
		overlappingGenes=new TreeSet<String>();
		for(String gene:gtf.mapGeneRange.keySet()) {
			Range range=gtf.mapGeneRange.get(gene);
			if(range.overlaps(pileup.from,pileup.to,pileup.seq)) {
				overlappingGenes.add(gene);
				numTranscripts+=gtf.mapGeneTranscripts.get(gene).size();
			}
		}
	}

	public void render(Pileup pileup, StringBuilder sb, int offsetY) {
		// TODO Auto-generated method stub
	
		
		
/*
		private void renderGTF(
				Pileup pileup, 
				StringBuilder sb, 
				GtfParser gtf, TreeSet<String> overlappingGenes, 
				double totalTrackWidth, double totalTrackHeight) {
			*/
			
//	renderGTF(sb, gtf, overlappingGenes, totalTrackWidth, totalTrackHeight);
	

	
		String colorFeatureLine="rgb(0,0,0)";

		//Write GTF annotation
		String textStyle="font: italic sans-serif; fill: red;";
		int curTrack=0;
		for(String gene:overlappingGenes) {
			for(String transcript:gtf.mapGeneTranscripts.get(gene)) {
				//System.out.println(transcript);
				//Should draw a straight line...

				double lineY=offsetY + (curTrack+0.5)*pileup.transcriptHeight;
				double featureY=lineY - pileup.featureHeight/2;
						
				//// Plot the line also suggesting the direction
				int minFrom=Integer.MAX_VALUE;
				int maxTo=0;
				for(Range r:gtf.mapTranscriptRanges.get(transcript)) {
					minFrom=Math.min(minFrom,r.from);
					maxTo=Math.max(maxTo,r.to);
				}
//				double minTransFrom=transformPos(minFrom);
	//			double maxTransTo=transformPos(maxTo);
				double minTransFrom=Math.max(pileup.labelsWidth,pileup.transformPosX(minFrom));
				double maxTransTo=Math.min(pileup.totalWidth,pileup.transformPosX(maxTo));

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
					double rFrom=pileup.transformPosX(r.from);
					double rTo=pileup.transformPosX(r.to);

					String exonColor="rgb(0,0,0)";
					if(r.featureType.equals(Range.FEATURE_3UTR))
						exonColor="rgb(0,255,0)";
					if(r.featureType.equals(Range.FEATURE_5UTR))
						exonColor="rgb(0,0,255)";
					
					String exonStyle="\"fill:"+exonColor+";stroke-width:3;stroke:none\"";
					
					sb.append("<rect x=\""+rFrom+"\" y=\""+featureY+"\" "
							+ "width=\""+(rTo-rFrom)+"\" "
							+ "height=\""+pileup.featureHeight+"\" style="+exonStyle+"/>");   //   mask=\"url(#trackmask)\"
					
					
					//TODO style should depend on exon, utr
				}

				//TODO: add a line, and direction
				
			//Add text: Transcript ID
			double textXFrom=5;
			double textY=lineY + pileup.textHeight/2;
			sb.append("<text x=\""+textXFrom+"\" y=\""+textY+"\" style=\""+textStyle+"\"  font-size=\""+pileup.textHeight+"px\" >"+transcript+"</text>");

			curTrack++;
			}
		}
	}

	@Override
	protected double getHeight(Pileup pileup) {
		return numTranscripts*pileup.transcriptHeight;
	}
	
	
	
}
