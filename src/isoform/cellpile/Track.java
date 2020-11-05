package isoform.cellpile;

public abstract class Track {

	public String trackName="";

	protected abstract double getHeight(Pileup pileup);

	protected abstract void render(Pileup pileup, StringBuilder sb, int offsetY);
	
	
}
