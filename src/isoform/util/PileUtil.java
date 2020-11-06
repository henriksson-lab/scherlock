package isoform.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.zip.GZIPInputStream;

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

	/**
	 * Clamp a value to a range
	 */
	public static int clamp(int x, int min, int max) {
		if(x<min)
			return(min);
		else if(x>max)
			return(max);
		else
			return(x);
	}

	
	/**
	 * Read list of zip file e.g. barcode list
	 */
	public static ArrayList<String> readBarcodeZipList(File fBarcodes) throws IOException {
		ArrayList<String> list=new ArrayList<String>();
		BufferedReader brFeatures=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fBarcodes))));
		String line;
		while((line=brFeatures.readLine())!=null) {
			list.add(line);
		}
		brFeatures.close();
		return list;
	}
	
}
