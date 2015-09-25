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
 * FastaIndexFileReader.java
 * Dec 5, 2014
 * 11:53:43 AM
 * yaping    lyping1986@gmail.com
 */
package edu.mit.compbio.qrf.java.utils;



import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.unc.genomics.Interval;

/**
 *
 */
public class FastaIndexFileReader {

	private Map<String,Integer> sequenceEntries = new HashMap<String,Integer>();
	private long genomeSize = 0;
	/**
	 * 
	 */
	public FastaIndexFileReader(File genomeIndex) {
		try {
			parseIndexFile(genomeIndex);
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	public boolean hasContig( String chr ) {
        return sequenceEntries.containsKey(chr);
    }
	
	public int size(){
		return sequenceEntries.size();
	}
	
	public long genomeSize(){
		return genomeSize;
	}
	
	public int getContigSize(String chr){
		if(!hasContig(chr)) throw new IllegalArgumentException("Can't find chromsome " + chr);
		return sequenceEntries.get(chr);
	}
	
	public Interval getContigAndLocation(long len){

		if(len > genomeSize) throw new IllegalArgumentException("Provided length " + len + " is bigger than genome length ");
		String tmp = null;
		for(String contig : sequenceEntries.keySet()){
			int size = sequenceEntries.get(contig);
			if(len - size < 0){
				return new Interval(contig, (int)len,(int)len);
			}
			len = len - size;
			tmp = contig;
		}
		return new Interval(tmp,(int) len,(int)len);
	}
	
	
	private void parseIndexFile(File genomeIndex) throws NumberFormatException, IOException{
		BufferedReader br = new BufferedReader(new FileReader(genomeIndex));
		String line;
		while( (line = br.readLine()) != null){
			String[] splitin = line.split("\t");
			sequenceEntries.put(splitin[0],Integer.parseInt(splitin[1]));
			genomeSize+=Long.parseLong(splitin[1]);
		}
		br.close();
		
	}
	

}
