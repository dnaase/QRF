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
 * VcfToMatrixQtlBinary.java
 * Dec 21, 2014
 * 10:30:31 PM
 * yaping    lyping1986@gmail.com
 */
package edu.mit.compbio.qrf.java.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import org.apache.log4j.Logger;

/**
 *
 */
public class VcfToMatrixQtlBinary extends HapmapToMatrixQtlBinary {

	protected static String USAGE = "VcfToMatrixQtlBinary [opts] input.hapmap.txt";

	protected static Logger log = Logger.getLogger(VcfToMatrixQtlBinary.class);	

	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		
		VcfToMatrixQtlBinary vmq = new VcfToMatrixQtlBinary();
		vmq.doMain(args);
	}	
	
	public VcfToMatrixQtlBinary() {
		super();
	}
	
	
	protected HashMap<String, Character> processLine(String inputHapmap,HashMap<String, Character> majorAlleleHash) throws IOException{
		String prefix = inputHapmap.replaceAll(".\\w+$", "");
		PrintWriter binaryMatrixWriter = new PrintWriter(prefix.concat(".matrixQtlBinary.txt"));
		
		
		HashMap<String, Character> majorAlleleOutputHash = new HashMap<String, Character>();
		
		log.info("Process input VCF file ...");
		BufferedReader br = new BufferedReader(new FileReader(inputHapmap));
		String line;
		
		int startCol = 10;
		while( (line = br.readLine()) != null){
			
			if(line.startsWith("##"))
				continue;
			String[] splitin = line.split("\t");
			
			if(splitin[0].startsWith("#CHROM")){
				binaryMatrixWriter.print("rs#");
				for(int i = startCol-1; i < splitin.length; i++){
					binaryMatrixWriter.print("\t" + splitin[i]);
				}
			}else{
				if(!splitin[6].equalsIgnoreCase("PASS"))
					continue;
				
				//if(splitin[2].equalsIgnoreCase(".")){
					binaryMatrixWriter.print("rs:" + splitin[0] + "-" + splitin[1]);
				//}else{
				//	binaryMatrixWriter.print(splitin[2]);
				//}
				
				
				//change vcf style into hapmap style
				HashMap<Character, Character> InfoToCharacter = new HashMap<Character, Character>();
				InfoToCharacter.put('.', 'N');
				InfoToCharacter.put('0', splitin[3].charAt(0));
				Integer s=1;
				for(String a : splitin[4].split(",")){
					InfoToCharacter.put(s.toString().charAt(0), a.charAt(0));
					s++;
				}
				//System.err.println(InfoToCharacter.get('.'));
				//System.err.println(InfoToCharacter.get('0'));
				//System.err.println(InfoToCharacter.get('1'));
				for(int i = startCol-1; i < splitin.length; i++){
					String[] infos = splitin[i].split(":");
					String[] genotypes = infos[0].split("/");
					splitin[i] = "";
					for(String genotype : genotypes){
						splitin[i] += InfoToCharacter.get(genotype.charAt(0));
						//System.err.println(splitin[1] + "\t" + genotype.charAt(0) + "\t" + InfoToCharacter.get(genotype.charAt(0)));
					}
				}
				//use hapmap style to calculate
				char majorAllele;
				if(majorAlleleHash != null && majorAlleleHash.containsKey(splitin[2])){
					majorAllele = majorAlleleHash.get(splitin[2]);
					majorAlleleOutputHash.put(splitin[2], majorAllele);
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
					majorAllele = majorAlleleSummary(splitin, startCol);
					majorAlleleOutputHash.put(splitin[0], majorAllele);
					//System.err.println(splitin[1] + "\t" + majorAllele);
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
	

	protected char majorAlleleSummary(String[] splitin, int startCol){
		HashMap<Character, Integer> alleleSummary = new HashMap<Character, Integer>();
		//String[] splitin = line.split("\t");
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

}
