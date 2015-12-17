/**
 * AnnotateBedWithHicBed.java
 * Nov 30, 2015
 * 11:31:23 AM
 * yaping    lyping1986@gmail.com
 */
package edu.mit.compbio.qrf.java.utils;

import htsjdk.tribble.Tribble;
import htsjdk.tribble.exception.UnsortedFileException;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;













import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.unc.genomics.BedGraphEntry;
import edu.unc.genomics.IntervalException;
import edu.unc.genomics.io.BedGraphFileReader;

/**
 *
 */
public class AnnotateBedWithHicBed {

	
	@Option(name="-sample",usage="specify the sample number in the VCF file if contain multiple sample genotype. number started at 1. default: 1")
	public int sample = 1;


	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "AnnotateBedWithHicBed [opts] input_hic.bed input.bed output.addHic.bed";
	
	private static Logger log = Logger.getLogger(AnnotateBedWithHicBed.class);

	private static long startTime = -1;

	private PrintWriter writer = null; 

	private final static int PURGE_INTERVAL = 100000;

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		AnnotateBedWithHicBed abwh = new AnnotateBedWithHicBed();
		BasicConfigurator.configure();
		abwh.doMain(args);

	}
	
	public void doMain(String[] args)
			throws Exception {

					CmdLineParser parser = new CmdLineParser(this);
					//parser.setUsageWidth(80);
					try
					{
						if(help || args.length < 3) throw new CmdLineException(USAGE);
						parser.parseArgument(args);
						
					
					}
					catch (CmdLineException e)
					{
						System.err.println(e.getMessage());
						// print the list of available options
						parser.printUsage(System.err);
						System.err.println();
						return;
					}

					String hicFile = arguments.get(0);
					String inputBedFile = arguments.get(1);
					String outputBedFile = arguments.get(2);
					BedGraphFileReader hicReader = initiate(hicFile, outputBedFile);
					
					//read input bed file, for each row, extract related HiC signal, write orignal line + HiC signal to the output bed file
					log.info("Parsing input bed file ...");
					BufferedReader br = new BufferedReader(new FileReader(inputBedFile));
					String line;
					long lineNum=0;

					while( (line = br.readLine()) != null){
						if(line.startsWith("#"))
							continue;
						String[] splitin = line.split("\t");
						String chr = splitin[0];
						int start = Integer.parseInt(splitin[1]);
						int end = Integer.parseInt(splitin[2]);
						if(end < start)
							throw new IntervalException("End is smaller than start coordinate. start: " + start + "\tend: " + end);
						
						Iterator<BedGraphEntry> overlapList = hicReader.query(chr, start, end);
						while(overlapList.hasNext()){
							
						}
						
						lineNum++;
						if(lineNum % PURGE_INTERVAL == 0){
							log.info("Processing line: " + lineNum);
						}
					}
					br.close();
					
					finish();
	}
	
	private BedGraphFileReader initiate(String hicFile, String outputFile) throws IOException{
		startTime = System.currentTimeMillis();
		
		
		writer = new PrintWriter(outputFile);
		
		
		 return new BedGraphFileReader(new File(hicFile).toPath());
	}

	private void finish() throws IOException{
		writer.close();
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;
		

		log.info("AnnotateBedWithHicBed's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}

}
