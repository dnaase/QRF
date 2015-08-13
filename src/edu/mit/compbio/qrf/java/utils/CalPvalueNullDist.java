/**
 * CalPvalueNullDist.java
 * Aug 13, 2015
 * 2:43:55 PM
 * yaping    lyping1986@gmail.com
 */
package edu.mit.compbio.qrf.java.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;


/**
 *
 */
public class CalPvalueNullDist {

	/**
	 * @param args
	 */
	
	@Option(name="-permutation",usage="number of permutation, default: 1000")
	public int permutation = 1000;
	
	@Option(name="-dataColRealPair",usage="the column number which contain values in the real pairs bed file, 1-based coordinate, default: 5")
	public int dataColRealPair = 5;
	
	@Option(name="-dataColNullPair",usage="the column number which contain values in the null distributed pairs bed file, 1-based coordinate, default: 5")
	public int dataColNullPair = 5;

	@Option(name="-genomicRangeStart",usage="the lower end of genomic range to calculate statistics test, default: 10000")
	public int genomicRangeStart = 10000;

	@Option(name="-genomicRangeEnd",usage="the upper end of genomic range to calculate statistics test, default: 100000")
	public int genomicRangeEnd = 100000;

	
	final private static String USAGE = "CalPvalueNullDist [opts] realPairs.bed nullPairs.bed";

	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	private PrintWriter writer = null; 
	private static final Logger log = Logger.getLogger(CalPvalueNullDist.class);
	
	private static long startTime = -1;
	private MersenneTwister generator = null;

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		CalPvalueNullDist cpnd = new CalPvalueNullDist();
		BasicConfigurator.configure();
		cpnd.doMain(args);

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
			if (arguments.size() < 2) throw new CmdLineException(USAGE);
		
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

		String realPairsBedFile = arguments.get(0);
		String nullPairsBedFile = arguments.get(1);
		
		initiate(realPairsBedFile);
		
		
		log.info("Parsing realPairs bed file ...");
		BufferedReader br = new BufferedReader(new FileReader(realPairsBedFile));
		String line;
		long lineNum=0;
		ArrayList<Double> realPairsValueList = new ArrayList<Double>();
		while( (line = br.readLine()) != null){
			if(line.startsWith("#"))
				continue;
			String[] splitin = line.split("\t");
			int start = Integer.parseInt(splitin[1]);
			int end = Integer.parseInt(splitin[2]);
			if((end - start) >= genomicRangeStart && (end - start) < genomicRangeEnd){
				realPairsValueList.add(Double.parseDouble(splitin[dataColRealPair]));
			}
			
			lineNum++;
			if(lineNum % 100000 == 0){
				log.info("Processing line: " + lineNum);
			}
		}
		br.close();
		double realPairsValueMedian = getMedian(realPairsValueList);
		log.info("realPairs median: " + realPairsValueMedian);
		
		log.info("Parsing nullPairs bed file ...");
		BufferedReader brNull = new BufferedReader(new FileReader(nullPairsBedFile));
		lineNum=0;
		ArrayList<Double> nullPairsValueList = new ArrayList<Double>();
		while( (line = brNull.readLine()) != null){
			if(line.startsWith("#"))
				continue;
			String[] splitin = line.split("\t");
			int start = Integer.parseInt(splitin[1]);
			int end = Integer.parseInt(splitin[2]);
			if((end - start) >= genomicRangeStart && (end - start) < genomicRangeEnd){
				if(!splitin[dataColNullPair].equalsIgnoreCase("NA") && !splitin[dataColNullPair].equalsIgnoreCase("NaN")){
					nullPairsValueList.add(Double.parseDouble(splitin[dataColNullPair]));
				}
				
			}
			lineNum++;
			if(lineNum % 100000 == 0){
				log.info("Processing line: " + lineNum);
			}
		}
		brNull.close();

		log.info("Sampling equal number of pairs from null distribution ...");
		lineNum=0;
		ArrayList<Double> nullPairsValueMedianList = new ArrayList<Double>();
		while(lineNum < permutation){
			ArrayList<Double> subSamplingNullPairsValueList = getSubSamplingList(nullPairsValueList, realPairsValueList.size());
			//log.info(subSamplingNullPairsValueList.get(0) + "\t" + subSamplingNullPairsValueList.get(1) + "\t" +subSamplingNullPairsValueList.get(2) + "\t" +subSamplingNullPairsValueList.get(3));
			nullPairsValueMedianList.add(getMedian(subSamplingNullPairsValueList));
			//log.info(getMedian(subSamplingNullPairsValueList));
			lineNum++;
			if(lineNum % 10 == 0){
				log.info("Permutation: " + lineNum);
			}
		}
		
		log.info("Calculate p value ...");
		int pos = findPosition(nullPairsValueMedianList, realPairsValueMedian, 0, nullPairsValueMedianList.size()-1);
		double permutatedPvalue = (double)(pos+1)/(double)nullPairsValueMedianList.size();
		log.info("p value is: " + permutatedPvalue);
		writer.println("realPairs median position is " + pos);
		writer.println("realPairs median is " + realPairsValueMedian);
		//log.info(nullPairsValueMedianList.get(0) + "\t" + nullPairsValueMedianList.get(1) + "\t" +nullPairsValueMedianList.get(2) + "\t" +nullPairsValueMedianList.get(3));
		writer.println("the median of null_pair_medians is " + getMedian(nullPairsValueMedianList));
		writer.println("p value is " + permutatedPvalue);
		finished();
	
	}
	
	private int findPosition(ArrayList<Double> permutatedProbList, double observedProb, int lower, int upper){
		if(lower == 0 && observedProb <= permutatedProbList.get(0)){
			return lower;
		}else if(upper == permutatedProbList.size()-1 && observedProb > permutatedProbList.get(permutatedProbList.size()-1)){
			return upper;
		}else if((upper-lower)==1 && (observedProb>=permutatedProbList.get(lower) && observedProb<=permutatedProbList.get(upper))){
			return upper;
		}else if(observedProb>=permutatedProbList.get(lower) && observedProb<=permutatedProbList.get(upper)){
			return findPosition(permutatedProbList, observedProb, lower+(upper-lower)/2, upper);
		}else if(observedProb < permutatedProbList.get(lower)){
			return findPosition(permutatedProbList, observedProb, Math.max(lower-(upper-lower),0), lower);
		}else{
			return findPosition(permutatedProbList, observedProb, upper, Math.max(upper+(upper-lower),permutatedProbList.size()-1));
		}
		
	}
	
	private double getMedian(ArrayList<Double> valueList){
		DescriptiveStatistics stats = new DescriptiveStatistics();
		for(double v : valueList){
			stats.addValue(v);
		}
		
		return stats.getPercentile(50);
		
	}
	
	private ArrayList<Double> getSubSamplingList(ArrayList<Double> valueList, int sampleSize){
		HashSet<Integer> randomIndexSet = new HashSet<Integer>();
		do{
			randomIndexSet.add(generator.nextInt(valueList.size()));
		}while(randomIndexSet.size()<sampleSize);
		
		ArrayList<Double> subSamplingList = new ArrayList<Double>();
		for(int i : randomIndexSet){
			subSamplingList.add(valueList.get(i));
		}
		
		
		return subSamplingList;
	}
	
	private void initiate(String outputFile) throws IOException{
		startTime = System.currentTimeMillis();
		String prefix = outputFile.replaceAll(".\\w+$", ".range-" + genomicRangeStart + "-" + genomicRangeEnd);
		dataColRealPair--;
		dataColNullPair--;
		writer = new PrintWriter(new File(prefix.concat(".permutationPvalue.txt")));
		generator = new MersenneTwister();
	}
	
	private void finished() throws IOException{

		writer.close();
		
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;


		log.info("CalPvalueNullDist's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");

	}

}
