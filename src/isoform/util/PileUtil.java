package isoform.util;

/**
 * Various utility functions
 * 
 * @author Johan Henriksson
 *
 */
public class PileUtil {

	/**
	 * Scales all the values to make max(list)=1
	 */
	public static void scaleMaxValue(double[] list, double scaleTo) {
		double theMax=PileUtil.maxForList(list);
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

	/**
	 * Return a string, repeated
	 */
	public static String[] getRepeatedString(String el, int n) {
		String[] list=new String[n];
		for(int i=0;i<n;i++)
			list[i] = el;
		return list;
	}

	public static int clamp(int x, int min, int max) {
		if(x<min)
			x=min;
		if(x>max)
			x=max;
		return(x);
	}

}
