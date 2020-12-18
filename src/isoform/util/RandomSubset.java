package isoform.util;
import java.io.File;
import java.io.IOException;
import java.util.Set;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * 
 * Subset reads from a SAM file
 * 
 * @author Johan Henriksson
 *
 */
public class RandomSubset {
	

	public static void subset(File fBAM, File fOut, Set<String> barcodes, double pKeepBC, double pKeepNonBC) throws IOException {
		//Open BAM in and output files
		final SamReader samreader = SamReaderFactory.makeDefault().open(fBAM);
		SAMFileHeader header = samreader.getFileHeader().clone();		
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		SAMFileWriter samwriter = factory.makeSAMOrBAMWriter(header, true,fOut);
		
		//To know progress, status variables
		int readRecords=0;
		int keptRecordsBC=0;
		int keptRecordsNonBC=0;

		//Loop through all SAM records
		for (final SAMRecord samRecord : samreader) {
			
			//Update user about progress
			readRecords++;
			if(readRecords%1000000 == 0){
				LogUtil.formatColumns(System.out, 25, 
						"Kept BC: "+keptRecordsBC,
						"Kept NonBC: "+keptRecordsNonBC,
						"Total: "+readRecords);
			}
				
			//Get BC for this read
			String bcCellCurrentCellBarcode=(String)samRecord.getAttribute("CB");
			
			//Check if this is a cell to count
			if(bcCellCurrentCellBarcode!=null && barcodes.contains(bcCellCurrentCellBarcode)) {
				if(Math.random()<pKeepBC) {
					keptRecordsBC++;
					samwriter.addAlignment(samRecord);
				}
			} else {
				if(Math.random()<pKeepNonBC) {
					keptRecordsNonBC++;
					samwriter.addAlignment(samRecord);
				}
			}
		}
			
		//Get out the last ones from memory. TODO remember to fix this in countBAM too
		samreader.close();
		samwriter.close();
		
		LogUtil.formatColumns(System.out, 25, 
				"Kept BC: "+keptRecordsBC,
				"Kept NonBC: "+keptRecordsNonBC,
				"Total: "+readRecords);
	}
	
	
	
	
	
}
