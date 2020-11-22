package isoform.trackRenderer;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.InputStream;
import java.util.Iterator;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.Interval;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixIndex;

// Note: A bit unclear if this really is a Tabix or if its a Tribble file format.
// 		 However, user not supposed to have to touch this, so doesnt matter that much.
public class TabixMaker {  
	

	public static String bgzipBedFile(String fileIn) throws Exception {

		String fileOut = fileIn + ".gz";
		InputStream is = new java.io.FileInputStream(fileIn);
		BufferedInputStream bis = new BufferedInputStream(is);
		BlockCompressedOutputStream bcos = new BlockCompressedOutputStream(fileOut);

		int bit = -1;
		while ((bit = bis.read()) != -1) {
			bcos.write(bit);
		}

		bcos.close();
		bis.close();
		is.close();
		
		return fileOut;
	}

	
	public static String indexBedFile(String bedFileName) throws Exception {
		final File bedFile = new File(bedFileName);
		final TabixIndex tabixIndex = IndexFactory.createTabixIndex(bedFile, new BEDDetailCodec(), null);
		tabixIndex.writeBasedOnFeatureFile(bedFile);
		final String indexFileName = Tribble.tabixIndexFile(bedFileName);
		return indexFileName;
	}
	
	
	// Test Tabix/Tribble file making
	public static void main(String[] args) throws Exception {
		
		System.out.println("main() starting...");

		// Make bgzipped file. Block gzipped files allows indexing and accessing portions of the data.
		// Thus, tabix wants this.
//		// Didnt work when I used file ending .bgzip, but .gz does.
//		String bedFileName = "C:\\Users\\anton\\java_projects\\isocounter\\data\\beds\\UP000005640_9606_lipid.bed";
		String bedFileName = "C:\\Users\\anton\\java_projects\\isocounter\\data\\beds\\merged_interesting.sorted.bed_details_compatible.chr_is_plain_number.bed";
		String bedBgzipFileName = TabixMaker.bgzipBedFile(bedFileName);
		
		// Make tabix index for bed file
//		String indexFileName = TabixMaker.indexBedFile(bedFileName);
//		System.out.println(indexFileName);
		String indexBgzipFileName = TabixMaker.indexBedFile(bedBgzipFileName);
//		System.out.println(indexBgzipFileName);

		// Query Tabix		
//		System.out.println(Tribble.indexFile(bedFileName));
//		System.out.println(Tribble.tabixIndexFile(bedFileName));
//		System.out.println(Tribble.indexFile(bedBgzipFileName));
//		System.out.println(Tribble.tabixIndexFile(bedBgzipFileName));
		
		final FeatureReader<BEDFeature> reader = AbstractFeatureReader.getFeatureReader(bedBgzipFileName, new BEDDetailCodec());

		// cd55 location
		Interval interval = new Interval("1", 207321376, 207360966);
		
		final Iterator<BEDFeature> readerIterator = 
				reader.query(interval.getContig(), interval.getStart(), interval.getEnd());
		while (readerIterator.hasNext()) {
			BEDFeature bedFeature = readerIterator.next();
//			System.out.println(bedFeature);
//			System.out.println(bedFeature.getColor());
			System.out.println(bedFeature.getDescription());
//			System.out.println(bedFeature.getExons());
			System.out.println(bedFeature.getLink());
			System.out.println(bedFeature.getName());
//			System.out.println(bedFeature.getScore());
//			System.out.println(bedFeature.getStrand());
//			System.out.println(bedFeature.getType());
			System.out.println(bedFeature.getStart());
			System.out.println(bedFeature.getEnd());
			System.out.println();
		}
		
		System.out.println("main() done!");
	}
}

