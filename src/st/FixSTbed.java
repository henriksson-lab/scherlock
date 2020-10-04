package st;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.StringTokenizer;

public class FixSTbed {
	public static void main(String[] args) throws IOException {
		
		File fIn=new File(args[0]);
		File fOut=new File(args[1]);
		
		BufferedReader br=new BufferedReader(new FileReader(fIn));
		PrintWriter os=new PrintWriter(new FileWriter(fOut));
		
		String line;
		while((line=br.readLine())!=null) {
			StringTokenizer stok=new StringTokenizer(line,"\t");
			String seq=stok.nextToken();
			int from=Integer.parseInt(stok.nextToken());
			int to=Integer.parseInt(stok.nextToken());
			os.println(seq+"\t"+Math.min(from, to)+"\t"+Math.max(from, to));
		}
		os.close();
		br.close();
	}

}
