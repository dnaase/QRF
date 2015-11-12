/**
 * HapmapToMatrixQtlBinary.java
 * Dec 20, 2014
 * 7:47:20 PM
 * yaping    lyping1986@gmail.com
 * 
 * convert Hapmap downloaded file (e.g. genotypes_chr1_CEU_r27_nr.b36_fwd.txt) into the input SNP file of matrixQTL
 * 0 will be major allele (if major allele frequency is 0.5, it will use reference allele), 1 to be the second allele.
 * if provide chimpanzze txt file, it will output one list of major allele have the same allele as chimpanzee, 
 * the other list: major allele do not agree with chimpazee. vague one will not be output
 */
package edu.mit.compbio.qrf.java.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

;

/**
 *
 */
public class HapmapToMatrixQtlBinary {

	@Option(name="-snp_id",usage="whether or not to output also SNP id map information in a seperate file. default:false")
	public boolean snp_id = false;
	@Option(name="-addChr",usage="whether or not to add chr as prefix in the contig number. default:false")
	public boolean addChr = false;
	@Option(name="-name",usage="column number of name field (rs ID) in snp138OrthoPt4Pa2Rm3.txt.1-based default:4")
	public int name = 4;
	@Option(name="-chimp_allele",usage="column number of chimpAllele field (e.g. chimpanzee's allele info) in snp138OrthoPt4Pa2Rm3.txt.1-based default:8")
	public int chimp_allele = 8;	
	@Option(name="-chimp_file",usage="File store chimpanzee allele information. snp138OrthoPt4Pa2Rm3.txt. default:null")
	public String chimp_file = null;	
	@Option(name="-dbsnp",usage="dbSNP.vcf file which contain major allele frequency. default:null")
	public String dbsnp = null;
	
	@Argument
	private List<String> arguments = new ArrayList<String>();
	

	protected static String USAGE = "HapmapToMatrixQtlBinary [opts] input.hapmap.txt";

	protected static Logger log = Logger.getLogger(HapmapToMatrixQtlBinary.class);	
	protected static long startTime = -1;
	protected static long lineNum=0;

	PrintWriter snpIdWriter = null;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		HapmapToMatrixQtlBinary hmq = new HapmapToMatrixQtlBinary();
		hmq.doMain(args);
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
						if (arguments.size() < 1) throw new CmdLineException(USAGE);
					
					}
					catch (CmdLineException e)
					{
						System.err.println(e.getMessage());
						// print the list of available options
						parser.printUsage(System.err);
						System.err.println();
						return;
					}
					//read input bed file, for each row,;
					String inputHapmap = arguments.get(0);
					initiate(inputHapmap);
					//get a hashmap of SNP name and major allele's character 
					HashMap<String, Character> majorAlleleHash = null;
					if(dbsnp != null){
						log.info("Input dbSNP major allele frequency hashmap");
						majorAlleleHash = getMajorAlleleHash(dbsnp);
					}else{
						log.info("Not provide input dbSNP.vcf, it will use input hashmap's own major allele");
					}
					//read into the
					log.info("Process input file ...");
					majorAlleleHash = processLine(inputHapmap, majorAlleleHash);
					
					//if chimp_file is not null, get a hashmap of SNP name and chimpanzee allele's character
					if(chimp_file != null){
						log.info("Get ancestal allele for each row in hashmap");
						getAncestalAllele(inputHapmap, chimp_file, majorAlleleHash);
					}
					
					finished();
	}
	
	protected HashMap<String, Character> getMajorAlleleHash(String dbsnpFile) throws IOException{
		HashMap<String, Character> majorAlleleHash =  new HashMap<String, Character>();
		BufferedReader br = new BufferedReader(new FileReader(dbsnpFile));
		String line;
		
		while( (line = br.readLine()) != null){
			if(line.startsWith("#"))
				continue;
			String[] splitin = line.split("\t");
			String keyName = splitin[2];
			String majorAllele;
			String[] altAlleles = splitin[4].split(",");
			String reg = "CAF=(\\d.\\d+)";
			for(int i = 0; i < altAlleles.length; i++){
				reg += ",(\\d.\\d+)";
			}
			Pattern pattern = Pattern.compile(reg);
			Matcher matcher = pattern.matcher(splitin[7]);
			if (matcher.find()){
				
				majorAllele = splitin[3];
				double max = Double.parseDouble(matcher.group(1));
				//System.err.println( splitin[4] + "\t" + splitin[7]);
				for(int i = 2, j = 0; i <= altAlleles.length+1;  i++, j++){
					
					double a = Double.parseDouble(matcher.group(i));
					if(a > max){
						majorAllele = altAlleles[j];
					}
				}
			}else{
				majorAllele = splitin[3];
				
			}
			majorAlleleHash.put(keyName, majorAllele.charAt(0));
			
		}
		br.close();
		
		return majorAlleleHash;
	}

	protected void getAncestalAllele(String inputHapmap, String chimpFile,HashMap<String, Character> majorAlleleHash) throws IOException{
		String prefix = inputHapmap.replaceAll(".\\w+$", "");
		
		PrintWriter ancestorListWriter = null;
		PrintWriter nonAncestorListWriter = null;

		ancestorListWriter = new PrintWriter(prefix.concat(".ancestorList.txt"));
		nonAncestorListWriter = new PrintWriter(prefix.concat(".nonAncestorList.txt"));
		
		BufferedReader br = new BufferedReader(new FileReader(chimpFile));
		String line;
		
		while( (line = br.readLine()) != null){
			if(line.startsWith("#") || line.startsWith("chrom"))
				continue;
			String[] splitin = line.split("\t");
			String keyName = splitin[name-1];
			char chimpAllele = splitin[chimp_allele-1].charAt(0);
			if(majorAlleleHash.containsKey(keyName)){
				char majorAllele = majorAlleleHash.get(keyName);
				if(majorAllele == chimpAllele){
					ancestorListWriter.println(keyName + "\t" + majorAllele + "\t" + chimpAllele);
				}else{
					nonAncestorListWriter.println(keyName + "\t" + majorAllele + "\t" + chimpAllele);
				}
			}
		}
		br.close();
		ancestorListWriter.close();
		nonAncestorListWriter.close();
		
	}
	
	//could be overwrite for the process of VCF file.
	protected HashMap<String, Character> processLine(String inputHapmap,HashMap<String, Character> majorAlleleHash) throws IOException{
		String prefix = inputHapmap.replaceAll(".\\w+$", "");
		PrintWriter binaryMatrixWriter = new PrintWriter(prefix.concat(".matrixQtlBinary.txt"));
		
		
		HashMap<String, Character> majorAlleleOutputHash = new HashMap<String, Character>();
		
		log.info("Process input file ...");
		BufferedReader br = new BufferedReader(new FileReader(inputHapmap));
		String line;
		
		int startCol = 12;
		while( (line = br.readLine()) != null){
			
			if(line.startsWith("#"))
				continue;
			String[] splitin = line.split(" ");
			binaryMatrixWriter.print(splitin[0]);
			if(splitin[0].startsWith("rs#")){
				for(int i = startCol-1; i < splitin.length; i++){
					binaryMatrixWriter.print("\t" + splitin[i]);
				}
				if(snp_id){
					snpIdWriter.println("SNP\tchrm_snp\tpos");
				}
			}else{
				if(snp_id){
					snpIdWriter.print(splitin[0]);
					if(addChr){
						snpIdWriter.println("\tchr" + splitin[2] + "\t" + splitin[3]);
					}else{
						snpIdWriter.println("\t" + splitin[2] + "\t" + splitin[3]);
					}
				}
				char majorAllele;
				if(majorAlleleHash != null && majorAlleleHash.containsKey(splitin[0])){
					majorAllele = majorAlleleHash.get(splitin[0]);
					majorAlleleOutputHash.put(splitin[0], majorAllele);
					for(int i = startCol-1; i < splitin.length; i++){
						char[] alleles = splitin[i].toCharArray();
						if(alleles[0] == 'N'){
							binaryMatrixWriter.print("\tNA");
						}else{
							Integer tmp = 0;
							for(char a : alleles){
								if(a != majorAllele)
									tmp++;
							}
							binaryMatrixWriter.print("\t" + tmp);
						}
					}
				}else{ //if there is no SNP information in dbSNP.vcf file, it will just use the summary statistics of the rows
					majorAllele = majorAlleleSummary(line, startCol);
					majorAlleleOutputHash.put(splitin[0], majorAllele);
					for(int i = startCol-1; i < splitin.length; i++){
						char[] alleles = splitin[i].toCharArray();
						if(alleles[0] == 'N'){
							binaryMatrixWriter.print("\tNA");
						}else{
							Integer tmp = 0;
							for(char a : alleles){
								if(a != majorAllele)
									tmp++;
							}
							binaryMatrixWriter.print("\t" + tmp);
						}
					}
				}
				
			}
			binaryMatrixWriter.println();
				
			lineNum++;
			if(lineNum % 5000 == 0){
				log.info("Processing line: " + lineNum);
			}
		}
		br.close();
		
		binaryMatrixWriter.close();
		
		return majorAlleleOutputHash;
	}
	
	protected char majorAlleleSummary(String line, int startCol){
		HashMap<Character, Integer> alleleSummary = new HashMap<Character, Integer>();
		String[] splitin = line.split(" ");
		for(int i = startCol-1; i < splitin.length; i++){
			char[] alleles = splitin[i].toCharArray();
			if(alleles[0] == 'N'){
				continue;
			}else{
				for(char a : alleles){
					int value = 1;
					if(alleleSummary.containsKey(a)){
						value = alleleSummary.get(a);
						value++;
					}
					alleleSummary.put(alleles[0], value);
				}
				
			}
		}
		int max = 0;
		char majorAllele = 'N';
		for(char a : alleleSummary.keySet()){
			int v = alleleSummary.get(a);
			if( v > max){
				max = v;
				majorAllele = a;
			}
		}
		return majorAllele;
	}
	
	private void initiate(String inputHapmap) throws IOException{
		startTime = System.currentTimeMillis();
		//System.out.println("HmmHunter started at : " + startTime);
		//System.out.println();
		if(snp_id){
			String prefix = inputHapmap.replaceAll(".\\w+$", "");
			snpIdWriter = new PrintWriter(prefix.concat(".snp_id_map.txt"));
		}

		
	}
	
	private void finished() throws IOException{
		if(snp_id){
			snpIdWriter.close();
		}
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;


		log.info( lineNum + " rows in total");
		log.info("HapmapToMatrixQtlBinary's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");

	}
	
	
}
