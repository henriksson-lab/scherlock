

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * 
 * @author Johan Henriksson
 *
 */
public class CountBam {

	
	//https://anaconda.org/bioconda/ucsc-gff3togenepred
	
	
	public static void main(String[] args) throws IOException {

		File fInput=new File("/home/mahogny/umeå/project/isoform/refgenome/Homo_sapiens.GRCh38.101.chr.gff3");
		if(args.length>0) {
			fInput=new File(args[0]);
		}
		System.out.println("doing "+fInput);

		File outdir=new File("./out");
		outdir.mkdirs();
		System.out.println("To: "+outdir);

		final SamReader reader = SamReaderFactory.makeDefault().open(fInput);

		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		HashMap<String, SAMFileWriter> libraryToWriter = new HashMap<String, SAMFileWriter>();

		int dups=0;
		int kept=0;
		int readRec=0;

		SAMRecord previousRecord=null;
		String bcCellPreviousU=null;
		String bcCellPreviousC=null;
		
		
		for (final SAMRecord samRecord : reader) {
			readRec++;
			if(readRec%10000 == 0){
				System.out.println("Now done "+readRec);
			}

			if(previousRecord!=null) {
				String bcCellCurrentU=(String)samRecord.getAttribute("UB");
				String bcCellCurrentC=(String)samRecord.getAttribute("CB");


				if(bcCellCurrentU!=null && bcCellCurrentC!=null) {

					if(bcCellCurrentU.equals(bcCellPreviousU) && bcCellCurrentC.equals(bcCellPreviousC)) {
						dups++;
					} else {
						//-F INT   only include reads with none of the FLAGS in INT present [0]
						//read unmapped (0x4)
						//not primary alignment (0x100)
						//read fails platform/vendor quality checks (0x200)
						//read is PCR or optical duplicate (0x400)
						//supplementary alignment (0x800)
						//in total: 0xF04

						int flags=samRecord.getFlags();
						if((flags & 0xF04) == 0) {
							//Only use as last if it had a barcode, and was of enough quality
							previousRecord=samRecord;

							String chrName=samRecord.getReferenceName();
							SAMFileWriter fw=libraryToWriter.get(chrName);
							if(fw==null) {
								SAMFileHeader header = reader.getFileHeader().clone();
								File fOut=new File(outdir,chrName+".bam");
								libraryToWriter.put(chrName, fw=factory.makeSAMOrBAMWriter(header, true,fOut));
							}

							fw.addAlignment(samRecord);
							kept++;
							
							bcCellPreviousU=bcCellCurrentU;
							bcCellPreviousC=bcCellCurrentC;

						}

					}
				}

			} else {
				previousRecord=samRecord;
			}


		}
		System.out.println("Kept: "+kept);
		System.out.println("Duplicates: "+dups);


		for(SAMFileWriter w:libraryToWriter.values()) {
			w.close();
		}
		reader.close();

	}

}

