package isoform.util;

/**
 * Various math utility functions
 * 
 * @author Johan Henriksson
 *
 */
public class PileMathUtil {

	/**
	 * Scales all the values to make max(list)=1
	 */
	public static void scaleMaxValue(double[] list, double scaleTo) {
		double theMax=PileMathUtil.maxForList(list);
		for(int i=0;i<list.length;i++) {
			list[i]*=scaleTo/theMax;
		}		
	}

	/**
	 * Get the maximum value in a list (double)
	 */
	public static double maxForList(double[] list) {
		double val=Double.MIN_VALUE;
		for(double x:list)
			val=Math.max(x,val);
		return val;
	}

	/**
	 * Get the maximum value in a list (int)
	 */
	public static int maxForList(int[] list) {
		int val=Integer.MIN_VALUE;
		for(int x:list)
			val=Math.max(x,val);
		return val;
	}

}
