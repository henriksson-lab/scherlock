package isoform.util;

import java.io.PrintStream;

/**
 * 
 * 
 * @author Johan Henriksson
 *
 */
public class LogUtil {
	
	
	public static void formatColumns(PrintStream os, int width, String... str) {
		
		//StringBuilder sb=new StringBuilder();
		for(String s:str) {
			os.print(s);
			for(int i=s.length();i<width;i++) {
				os.print(' ');
			}
		}
		os.println();
	}
	
	public static void main(String[] args) {
		formatColumns(System.out, 10,"foo","bar");
	}
	
	

}
