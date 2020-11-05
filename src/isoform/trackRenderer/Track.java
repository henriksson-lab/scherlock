package isoform.trackRenderer;

/**
 * 
 * A track that can be rendered
 * 
 * @author Johan Henriksson
 *
 */
public abstract class Track {

	public String trackName="";

	protected abstract double getHeight(TrackRenderer renderer);
	protected abstract void render(TrackRenderer renderer, StringBuilder sb, int offsetY);
	protected abstract void allocateSize(TrackRenderer renderer);

}
