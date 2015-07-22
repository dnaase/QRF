/*
 * The MIT License (MIT)
 * Copyright (c) 2015 dnaase <Yaping Liu: lyping1986@gmail.com>

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * QRF.java
 * Dec 9, 2014
 * 5:34:04 PM
 * yaping    lyping1986@gmail.com
 */
package edu.mit.compbio.qrf.java;



import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;

import net.sf.javaml.classification.Classifier;
import net.sf.javaml.classification.evaluation.CrossValidation;
import net.sf.javaml.classification.evaluation.PerformanceMeasure;
import net.sf.javaml.classification.tree.RandomForest;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.Instance;
import net.sf.javaml.tools.data.FileHandler;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;



public class QRF {


	/**
	 * @param args
	 */

	@Option(name="-trainOnly",usage="only  in the train mode, default: false")
	public boolean trainOnly = false;
	
	@Option(name="-ignoreCV",usage="ignore the cross validation step, default: false")
	public boolean ignoreCV = false;

//	@Option(name="-trainFile",usage="the train data set file name to be load, default: null")
//	public String trainFile = null;
	
	@Option(name="-inputFile",usage="the test data set file name to be load, default: null")
	public String testFile = null;
	
	@Option(name="-outputFile",usage="after apply model, the model predicted value on test dataset, default: null")
	public String testOutput = null;

	@Option(name="-seed",usage="seed for the randomization, default: 12345")
	public Integer seed = 12345;
	
	@Option(name="-kFold",usage="number of cross validation during training step, default: 10")
	public int kFold = 10;

	@Option(name="-numTrees",usage="number of tree for random forest model, default: 1000")
	public Integer numTrees = 1000;
	
//	@Option(name="-maxDepth",usage="maximum number of tree depth for random forest model, default: 4")
//	public Integer maxDepth = 4;
	
//	@Option(name="-maxBins",usage="maximum number of bins used for splitting features at random forest model, default: 100")
//	public Integer maxBins = 100;

	@Option(name="-classIndex",usage="which column is the class to be identified, default: 5")
	public int classIndex = 4;
	
	@Option(name="-permutation",usage="enable permutation test to calculate p value, default: not enabled")
	public boolean permutation = false;
	
	@Option(name="-permutationClass",usage="the class name used to do permutation,e.g: predict of meqtl/nomeqtl, use meqtl's probability to do permutation test, default: meqtl")
	public String permutationClass = "meqtl";

	@Option(name="-permutationNum",usage="permutation number. When the number of total line in testOutput is small (<1M lines), it will be better to make thie number larger. default: 1")
	public int permutationNum = 1;


	@Option(name="-sep",usage="seperate character to split each column, default: \t")
	public String sep = "\t";

	@Option(name="-h",usage="show option information")
	public boolean help = false;
	
	final private static String USAGE = "QRF [opts] trainFile.txt";

	@Argument
	private List<String> arguments = new ArrayList<String>();

	
	private static Logger log = Logger.getLogger(QRF.class);
	private PrintWriter writer = null; 
	private static long startTime = -1;

	
	private Dataset data = null;
	

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		QRF qrf = new QRF();
		BasicConfigurator.configure();
		qrf.doMain(args);
	}

	public void doMain(String[] args)
			throws Exception {

					CmdLineParser parser = new CmdLineParser(this);
					parser.setUsageWidth(80);
					try
					{
						if(help || args.length < 1) throw new CmdLineException(USAGE);
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
					String trainFile = args[0];
					//read input bed file, for each row,
					//String trainModelFile = arguments.get(0);
					initiate(trainFile, testOutput);
					
					
					Classifier model = null;
					if(trainFile != null){
						model = trainModel();
						//saveModel(model, trainModelFile);
					}else{
						throw new FileNotFoundException("no input train file");//model = loadModel(trainModelFile);
					}
					
					if(!trainOnly && testFile != null){
						testModel(model, testFile);
					}
					
					finish(testOutput);
	}
	
	private void initiate(String trainFile, String outputFile) throws IOException{
		startTime = System.currentTimeMillis();
		 data = FileHandler.loadDataset(new File(trainFile), classIndex-1, sep);
		if(outputFile != null){
			writer = new PrintWriter(new File(outputFile));
		}
		
	}

	private void finish(String outputFile) throws IOException{
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;
		
		
		log.info("testRandomForestClassifier's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}

	private Classifier trainModel(){
		log.info("Training model ...");
		//System.err.println(data.noAttributes());
		//Classifier rf = new RandomForest(numTrees);
		Classifier rf = new RandomForest(numTrees, false, 1, new Random(seed));
		
		if(!ignoreCV){
			//Construct new cross validation instance with the RandomForest classifier 
			CrossValidation cv = new CrossValidation(rf);
			//Perform cross-validation on the data set
			Map<Object, PerformanceMeasure> p = cv.crossValidation(data, kFold, new Random(seed));
			for(Object key : p.keySet()){
				PerformanceMeasure val = p.get(key);
				log.info(key + "\tError rate: " + val.getErrorRate() + "\tPrecision:" + val.getPrecision() + "\tRecall: " + val.getRecall() + "\tAUC: " + val.getBCR());
				log.info(val);
			}
			
		}
		

		//train on the whole data set		
		rf.buildClassifier(data);
		
		return rf;

	}
	
	private void testModel(Classifier model, String testFile) throws IOException{
		log.info("Testing model ...");
		 Dataset dataForClassification = FileHandler.loadDataset(new File(testFile), classIndex-1, sep);
		 
		 long num=0;
		 ArrayList<Double> probabilityList = new ArrayList<Double>();
		 /* Classify all instances and check with the correct class values */
		 for (Instance inst : dataForClassification) {
		     
			 Object predictedClassValue = model.classify(inst);
		     
		     for(Integer key : inst.keySet()){
		    	 writer.print(inst.get(key) + "\t");
		     }
		     writer.print(inst.classValue() + "\t");
		     writer.print(predictedClassValue.toString());
		     Map<Object, Double> probMap = model.classDistribution(inst);
		     if(permutation){
		    	 probabilityList.add(probMap.get(permutationClass));
		     }
		     writer.println();
		     num++;
		     if(num % 100000 == 0){
		    	 log.info("Line " + num + " ...");
		    	 writer.flush();
		     }
		 }
		 
		 if(testOutput != null){
				writer.close();
		}
		 
		if(permutation){
			permutationTest(model, dataForClassification, probabilityList);
		}
		
	}
	
	private void permutationTest(Classifier model,Dataset dataForClassification, ArrayList<Double> probabilityList) throws IOException{
		log.info("Permutation ...");
		HashSet<Double> permutatedProbSet = new HashSet<Double>();
		
		for(int j=1; j<=permutationNum; j++){
			Collections.shuffle(dataForClassification);
			long num=0;
			for (Instance inst : dataForClassification) {
				double tmp = model.classDistribution(inst).get(permutationClass);
				permutatedProbSet.add(tmp);
				num++;
			     if(num % 100000 == 0){
			    	 log.info("Line " + num + " ...");
			    	 writer.flush();
			     }
			}
			 log.info("Permutation in time  " + j + " ...");
		}
		
		
		log.info("Calculate permutate p value ...");
		//calculate permutate p value
		ArrayList<Double> permutatedProbList = new ArrayList<Double>(permutatedProbSet);
		Collections.sort(permutatedProbList);
		
		HashMap<Integer, String> pvalueAndLabel = new HashMap<Integer, String>();
		int num = 0;
		for(double observedProb : probabilityList){
			int pos = findPosition(permutatedProbList, observedProb, 0, permutatedProbList.size()-1);
			double permutatedPvalue = (double)(permutatedProbList.size()-pos)/(double)permutatedProbList.size();
			pvalueAndLabel.put(num, String.valueOf(permutatedPvalue));
			num++;
		}
		
		log.info("Move files ...");
		
		String tmpFile = testOutput + ".permuatetdTmp.txt";
		BufferedReader br = new BufferedReader(new FileReader(testOutput));
		PrintWriter tmpWriter = new PrintWriter(new File(tmpFile));
		
		String line;
		int i = 0;
		while( (line = br.readLine()) != null){
			if(line.startsWith("#")){
				tmpWriter.println(line + "\tpermutatedP");
				continue;
			}
			
			tmpWriter.println(line + "\t" + String.format("%.4f",Double.parseDouble(pvalueAndLabel.get(i))) );
			i++;
		}
		br.close();
		tmpWriter.close();
		Files.move(new File(tmpFile).toPath(), new File(testOutput).toPath(), StandardCopyOption.REPLACE_EXISTING);
		
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
}
