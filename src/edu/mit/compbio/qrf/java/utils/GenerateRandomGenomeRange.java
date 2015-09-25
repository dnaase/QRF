/**
 * GenerateRandomGenomeRange.java
 * Sep 25, 2015
 * 4:03:32 PM
 * yaping    lyping1986@gmail.com
 */
package edu.mit.compbio.qrf.java.utils;

import htsjdk.tribble.Tribble;
import htsjdk.tribble.exception.UnsortedFileException;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.readers.TabixReader.Iterator;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.unc.genomics.Interval;
import edu.unc.genomics.io.BedFileReader;
import edu.unc.genomics.io.BedGraphFileReader;
import edu.unc.genomics.io.IntervalFileReader;



/**
 *
 */
public class GenerateRandomGenomeRange {


	@Option(name="-iteration",usage="the maximum iteration allowed to search the best matched random reigon, default: 1000")
	public int iteration = 1000;
	@Option(name="-exclude",usage="the genomic region file to excluded in the random region, default: not enabled")
	public String exclude = null;
	@Option(name="-randomStart",usage="the random genomic region start , default: 1")
	public int start = 1;
	@Option(name="-randomEnd",usage="the random genomic region end , -randomStart 1 -randomEnd 1000 means the region length would be uniformly random in 1-1000 region. default: 1000")
	public int end = 1000;
	
	
	final private static String USAGE = "GenerateRandomGenomeRange [opts] interval_num output.bed hg19.fa.fai";

	@Argument
	private List<String> arguments = new ArrayList<String>();

	private PrintWriter bedWriter = null; 

	private File genomeIndex = null;
	private MersenneTwister generator = null;

	private RandomGenomicRegionGenerator randomGenomicGenerator = null;
	private static final Logger log = Logger.getLogger(GenerateRandomGenomeRange.class);
	
	private static long startTime = -1;
	private static long lineNum=0;
	
	private IntervalFileReader bedReader = null;

	
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		GenerateRandomGenomeRange grg = new GenerateRandomGenomeRange();
		BasicConfigurator.configure();
		grg.doMain(args);

	}
	
	public void doMain(String[] args)
	throws Exception {

			CmdLineParser parser = new CmdLineParser(this);
			// if you have a wider console, you could increase the value;
			// here 80 is also the default
			parser.setUsageWidth(80);
			try
			{
				parser.parseArgument(args);
				if (arguments.size() < 3) throw new CmdLineException(USAGE);
			
			}
			catch (CmdLineException e)
			{
				System.err.println(e.getMessage());
				// print the list of available options
				parser.printUsage(System.err);
				System.err.println();
				return;
			}
			//read input bed file, for each row,
			long randomNum = Long.parseLong(arguments.get(0));
			String outputBedFile = arguments.get(1);
			genomeIndex = new File(arguments.get(2));
			
			
			initiate(outputBedFile);
			int range = end-start;
			for(long i = 0; i < randomNum; i++){
				int randomBlockLen = generator.nextInt(range)+1;
				Interval randomRegion = null;
				int it = 0;
				if(exclude != null){
					do{
						if(it > iteration){
							randomBlockLen = generator.nextInt(range)+1;
						}
						randomRegion = randomGenomicGenerator.next(randomBlockLen);
						
						it++;
					}while(bedReader.query(randomRegion.getChr(), randomRegion.getStart(), randomRegion.getStop()).hasNext());
					
				}else{
					//System.err.println(randomBlockLen + "\t" + randomGenomicGenerator.randomRange);
					randomRegion = randomGenomicGenerator.next(randomBlockLen);
				}
				bedWriter.println(randomRegion.getChr() + "\t" + randomRegion.getStart() + "\t" + randomRegion.getStop() + "\t" + randomRegion.toOutput());
				lineNum++;
				if(lineNum % 5000 == 0){
					log.info("Processing line: " + lineNum);
				}
			}
			
			finished();
	}
	
	
	private void initiate(String outputFile) throws IOException{
		startTime = System.currentTimeMillis();
		//System.out.println("HmmHunter started at : " + startTime);
		//System.out.println();
		randomGenomicGenerator  = new RandomGenomicRegionGenerator(genomeIndex);
		generator = new MersenneTwister();
		
		bedWriter = new PrintWriter(new File(outputFile));
		
		if(exclude != null){
			
			bedReader = new BedFileReader(new File(exclude).toPath());
			
			//bedReader = IntervalFileReader.autodetect(new File(exclude).toPath());
		}
		
		
		 
		
		
	}
	
	private void finished() throws IOException{
		
		
		bedWriter.close();
		if(exclude != null){
			bedReader.close();
		}
		
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;


		log.info("GenerateRandomGenomeRange's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");

	}

}
