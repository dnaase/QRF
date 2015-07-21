# The MIT License (MIT)
# Copyright (c) 2015 dnaase <Yaping Liu: lyping1986@gmail.com>

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# matrixEqtlEngine.R
# Nov 26, 2014
# 10:09:44 AM
# 
# Author: yaping
###############################################################################

library(MatrixEQTL)


#####matrix eQTL without SNP/Expre location information
matrixEQtlAll<-function(snpInfo, exprInfo, covarInfo, qqplotFile="qqplot.matrixEQtlAll.pdf", outputFile="qtl.matrixEQtlAll.txt", pvalue=1){
	useModel = modelLINEAR;
	
	# Output file name
	output_file_name = outputFile;
	
	# Only associations significant at this level will be saved
	pvOutputThreshold = pvalue;
	
	# Error covariance matrix
	# Set to numeric() for identity.
	errorCovariance = numeric();
	
	## Load genotype data
	
	snps = SlicedData$new();
	snps$fileDelimiter = "\t";      # the TAB character
	snps$fileOmitCharacters = "NA"; # denote missing values;
	snps$fileSkipRows = 1;          # one row of column labels
	snps$fileSkipColumns = 1;       # one column of row labels
	snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
	snps$LoadFile(snpInfo);
	
	## Load gene expression data
	
	gene = SlicedData$new();
	gene$fileDelimiter = "\t";      # the TAB character
	gene$fileOmitCharacters = "NA"; # denote missing values;
	gene$fileSkipRows = 1;          # one row of column labels
	gene$fileSkipColumns = 1;       # one column of row labels
	gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
	gene$LoadFile(exprInfo);
	
	## Load covariates
	
	cvrt = SlicedData$new();
	cvrt$fileDelimiter = "\t";      # the TAB character
	cvrt$fileOmitCharacters = "NA"; # denote missing values;
	cvrt$fileSkipRows = 1;          # one row of column labels
	cvrt$fileSkipColumns = 1;       # one column of row labels
	if(length(covarInfo)>0) {
		cvrt$LoadFile(covarInfo);
	}
	
	## Run the analysis
	
	me = Matrix_eQTL_engine(
			snps = snps,
			gene = gene,
			cvrt = cvrt,
			output_file_name = output_file_name,
			pvOutputThreshold = pvOutputThreshold,
			useModel = useModel, 
			errorCovariance = errorCovariance, 
			verbose = TRUE,
			pvalue.hist = "qqplot",
			min.pv.by.genesnp = FALSE,
			noFDRsaveMemory = FALSE);
	
	#unlink(output_file_name);
	
	## Results:
	
	cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
	#cat('Detected eQTLs:', '\n');
	#show(me$all$eqtls)
	pdf(qqplotFile, paper="special", height=5, width=5)
	plot(me, pch=16, cex=0.7)
	dev.off()
	me
}


#####matrix eQTL with cis trans information
matrixEQtlCis<-function(snpInfo, exprInfo, covarInfo, qqplotFile="qqplot.matrixEQtlCis.pdf", snpLoc, exprLoc, outputFileCis="cis-qtl.matrixEQtlAll.txt", outputFileTrans="trans-qtl.matrixEQtlAll.txt", pvalueCis=1, pvalueTrans=1, cisDist=1e6){
	useModel = modelLINEAR;
	output_file_name_cis = outputFileCis;
	output_file_name_tra = outputFileTrans;
	
	# Only associations significant at this level will be saved
	pvOutputThreshold_cis = pvalueCis;
	pvOutputThreshold_tra = pvalueTrans;
	
	# Error covariance matrix
	# Set to numeric() for identity.
	errorCovariance = numeric();
	
	# Distance for local gene-SNP pairs
	cisDist = cisDist;
	
	## Load genotype data
	
	snps = SlicedData$new();
	snps$fileDelimiter = "\t";      # the TAB character
	snps$fileOmitCharacters = "NA"; # denote missing values;
	snps$fileSkipRows = 1;          # one row of column labels
	snps$fileSkipColumns = 1;       # one column of row labels
	snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
	snps$LoadFile(snpInfo);
	
	## Load gene expression data
	
	gene = SlicedData$new();
	gene$fileDelimiter = "\t";      # the TAB character
	gene$fileOmitCharacters = "NA"; # denote missing values;
	gene$fileSkipRows = 1;          # one row of column labels
	gene$fileSkipColumns = 1;       # one column of row labels
	gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
	gene$LoadFile(exprInfo);
	
	## Load covariates
	
	cvrt = SlicedData$new();
	cvrt$fileDelimiter = "\t";      # the TAB character
	cvrt$fileOmitCharacters = "NA"; # denote missing values;
	cvrt$fileSkipRows = 1;          # one row of column labels
	cvrt$fileSkipColumns = 1;       # one column of row labels
	if(length(covarInfo)>0) {
		cvrt$LoadFile(covarInfo);
	}
	
	## Run the analysis
	snpspos = read.table(snpLoc, header = TRUE, stringsAsFactors = FALSE);
	genepos = read.table(exprLoc, header = TRUE, stringsAsFactors = FALSE);
	
	me = Matrix_eQTL_main(
			snps = snps, 
			gene = gene, 
			cvrt = cvrt,
			output_file_name     = output_file_name_tra,
			pvOutputThreshold     = pvOutputThreshold_tra,
			useModel = useModel, 
			errorCovariance = errorCovariance, 
			verbose = TRUE, 
			output_file_name.cis = output_file_name_cis,
			pvOutputThreshold.cis = pvOutputThreshold_cis,
			snpspos = snpspos, 
			genepos = genepos,
			cisDist = cisDist,
			pvalue.hist = "qqplot",
			min.pv.by.genesnp = FALSE,
			noFDRsaveMemory = FALSE);
	
	#unlink(output_file_name_tra);
	#unlink(output_file_name_cis);
	
	## Results:
	
	cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
	#cat('Detected local eQTLs:', '\n');
	#show(me$cis$eqtls)
	
	#cat('Detected distant eQTLs:', '\n');
	#show(me$trans$eqtls)
	
	pdf(qqplotFile, paper="special", height=5, width=5)
	plot(me, pch=16, cex=0.7)
	dev.off()
	me
}

