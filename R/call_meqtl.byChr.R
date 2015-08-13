# call_meqtl.byChr.R
# Aug 6, 2015
# 12:29:04 PM
# 
# Author: yaping
###############################################################################

library(MatrixEQTL)

qrf_path="./"
wd="./"
covarInfo=NULL
chr="chr1"
for (e in commandArgs(TRUE)) {
	ta = strsplit(e,"=",fixed=TRUE)
	if(! is.na(ta[[1]][2])) {
		if(ta[[1]][1] == "qrf_path"){
			qrf_path<-ta[[1]][2] 
		}
		if(ta[[1]][1] == "wd"){
			wd<-ta[[1]][2] 
		}
		if(ta[[1]][1] == "prefix"){
			prefix<-ta[[1]][2]
		}
		if(ta[[1]][1] == "chr"){
			chr<-ta[[1]][2]
		}
		if(ta[[1]][1] == "snpInfo"){
			snpInfo<-ta[[1]][2]
		}
		if(ta[[1]][1] == "exprInfo"){
			exprInfo<-ta[[1]][2]
		}
		if(ta[[1]][1] == "snpLoc"){
			snpLoc<-ta[[1]][2]
		}
		if(ta[[1]][1] == "exprLoc"){
			exprLoc<-ta[[1]][2]
		}
		if(ta[[1]][1] == "covarInfo"){
			covarInfo<-ta[[1]][2]
		}
	}
}
source(paste(qrf_path,"/R/matrixEqtlEngine.R",sep=""))

setwd(wd)

snp<-read.table(snpInfo,sep="\t",header=T)
expr<-read.table(exprInfo,sep="\t",header=T)


snp_loc<-read.table(snpLoc,sep="\t",header=T)
expr_loc<-read.table(exprLoc,sep="\t",header=T)

snp_loc<-snp_loc[snp_loc[,2]%in%chr,]
expr_loc<-expr_loc[expr_loc[,2]%in%chr,]
colnames(snp_loc)[1]<-"SNP"
colnames(expr_loc)[1]<-"geneid"
snp_chr=paste(snpLoc,".",prefix,".",chr,".txt",sep="")
expr_chr=paste(exprLoc,".",prefix,".",chr,".txt",sep="")

write.table(snp_loc,snp_chr,sep="\t",quote =F, row.names=F,col.names =T)
write.table(expr_loc,expr_chr,sep="\t",quote =F, row.names=F,col.names =T)


snp_1<-snp[snp[,1] %in% snp_loc[,1],]
expr_1<-expr[expr[,1] %in% expr_loc[,1],]

colnames(snp_1)[1]<-"id"
colnames(expr_1)[1]<-"id"


snp_name=paste(snpInfo,".",prefix,".",chr,".txt",sep="")
expr_name=paste(exprInfo,".",prefix,".",chr,".txt",sep="")


write.table(snp_1,snp_name,sep="\t",quote =F, row.names=F,col.names =T)
write.table(expr_1,expr_name,sep="\t",quote =F, row.names=F,col.names =T)


qqplotFile=paste("qqplot.matrixEQtlCis.",prefix,".",chr,".pdf",sep="")
outputFileCis=paste("cis-qtl.matrixEQtlAll.",prefix,".",chr,".txt",sep="")
outputFileTrans=paste("trans-qtl.matrixEQtlAll.",prefix,".",chr,".txt",sep="")
meQTL.random<-matrixEQtlCis(snpInfo=snp_name, exprInfo=expr_name, covarInfo=covarInfo, qqplotFile=qqplotFile, snpLoc=snp_chr, exprLoc=expr_chr, outputFileCis=outputFileCis, outputFileTrans=outputFileTrans,pvalueCis=0.01, pvalueTrans=1e-300, pvalue.hist=FALSE)

file.remove(snp_name)
file.remove(expr_name)
file.remove(snp_chr)
file.remove(expr_chr)


