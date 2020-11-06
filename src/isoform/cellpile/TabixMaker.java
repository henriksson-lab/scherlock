package isoform.cellpile;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.InputStream;
import java.util.Iterator;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.Interval;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixIndex;


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
	
	
	public static void main(String[] args) throws Exception {

		// Make bgzipped file. Block gzipped files allows indexing and accessing portions of the data.
		// Thus, tabix wants this.
//		// Didnt work when I used file ending .bgzip, but .gz does.
//		String bedFileName = "C:\\Users\\anton\\java_projects\\isocounter\\data\\beds\\UP000005640_9606_lipid.bed";
		String bedFileName = "C:\\Users\\anton\\java_projects\\isocounter\\data\\beds\\merged_all.sorted.bed_details_compatible.bed";
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
		Interval interval = new Interval("chr1", 207321376, 207360966);
		
		final Iterator<BEDFeature> readerIterator = 
				reader.query(interval.getContig(), interval.getStart(), interval.getEnd());
		while (readerIterator.hasNext()) {
			BEDFeature bedFeature = readerIterator.next();
//			System.out.println(bedFeature);
//			System.out.println(bedFeature.getColor());
			System.out.println(bedFeature.getDescription());
//			System.out.println(bedFeature.getExons());
			System.out.println(bedFeature.getLink());
//			System.out.println(bedFeature.getName());
//			System.out.println(bedFeature.getScore());
//			System.out.println(bedFeature.getStrand());
//			System.out.println(bedFeature.getType());
//			System.out.println(bedFeature.getStart());
//			System.out.println(bedFeature.getEnd());
			System.out.println();
		}

	}

}


//cd55: "1",207321376,207360966



//// iterate over the query intervals and validate the query results
//try(final FeatureReader<BEDFeature> originalReader =
//          AbstractFeatureReader.getFeatureReader(inputBed.getAbsolutePath(), new BEDCodec());
//  final FeatureReader<BEDFeature> createdReader =
//          AbstractFeatureReader.getFeatureReader(tmpBed.getAbsolutePath(), new BEDCodec()))
//{
//  for (final Interval interval: queryIntervals) {
//      final Iterator<BEDFeature> originalIt = originalReader.query(interval.getContig(), interval.getStart(), interval.getEnd());
//      final Iterator<BEDFeature> createdIt = createdReader.query(interval.getContig(), interval.getStart(), interval.getEnd());
//      while(originalIt.hasNext()) {
//          Assert.assertTrue(createdIt.hasNext(), "some features not returned from query");
//          BEDFeature bedOrig = originalIt.next();
//          BEDFeature bedTmp = createdIt.next();
//          Assert.assertEquals(bedOrig.getContig(), bedTmp.getContig());
//          Assert.assertEquals(bedOrig.getStart(), bedTmp.getStart());
//          Assert.assertEquals(bedOrig.getEnd(), bedTmp.getEnd());
//      }
//  }
//}
