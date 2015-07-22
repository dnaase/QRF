##perl -e '@randomSizes=(575);foreach $random(@randomSizes){for($i=1;$i<=1;$i++){$prefix="matt_meqtl_random_$i";@chrs=(21,22);foreach $chr(@chrs){$chr="chr".$chr;`perl ~/qsubexc.long.pl "/broad/software/free/Linux/redhat_6_x86_64/pkgs/r_3.1.1-bioconductor-3.0/bin/R --no-save --no-restore --args wd=/home/unix/lyping/hptmp/AD_meqtl_tmp/partition_maskCG_all_chr/ snpInfo=SNP.tab exprInfo=DNAm.mask.tab  snpLoc=SNP.locs exprLoc=DNAm.mask.locs covarInfo=covs.tab prefix=$prefix chr=$chr sampleSize=$random < call_meqtl_random.byChr.R" ${prefix}_size_${chr}_$random 10 1`;}}}'

library(MatrixEQTL)

qrf_path="./"
wd="./"
covarInfo=NULL
chr="chr1"
sampleSize=NULL
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
		if(ta[[1]][1] == "sampleSize"){
			sampleSize<-as.numeric(ta[[1]][2])
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
if(!is.null(covarInfo)){
	cov<-read.table(covarInfo,sep="\t",header=T)
}

if(is.null(sampleSize)){
	sampleSize<-dim(snp)[2]-1
}


random_1_names<-sample(colnames(snp)[2:length(snp[1,])],sampleSize)


snp_loc<-read.table(snpLoc,sep="\t",header=T)
expr_loc<-read.table(exprLoc,sep="\t",header=T)

snp_loc<-snp_loc[snp_loc[,2]%in%chr,]
expr_loc<-expr_loc[expr_loc[,2]%in%chr,]
colnames(snp_loc)[1]<-"SNP"
colnames(expr_loc)[1]<-"geneid"
snp_chr=paste(snpLoc,".SampleSize-",sampleSize,".",prefix,".",chr,".txt",sep="")
expr_chr=paste(exprLoc,".SampleSize-",sampleSize,".",prefix,".",chr,".txt",sep="")

write.table(snp_loc,snp_chr,sep="\t",quote =F, row.names=F,col.names =T)
write.table(expr_loc,expr_chr,sep="\t",quote =F, row.names=F,col.names =T)



snp_1<-snp[,colnames(snp) %in% random_1_names]
expr_1<-expr[,colnames(expr) %in% random_1_names]
if(!is.null(covarInfo)){
	cov_1<-cov[,colnames(cov) %in% random_1_names]
	cov_1<-cbind(cov[,1],cov_1)
	colnames(cov_1)[1]<-"id"
	cov_random=paste(covarInfo,".SampleSize-",sampleSize,".",prefix,".",chr,".txt",sep="")
	write.table(cov_1,cov_random,sep="\t",quote =F, row.names=F,col.names =T)
}else{
	cov_random=NULL
}

snp_1<-cbind(snp[,1],snp_1)
expr_1<-cbind(expr[,1],expr_1)


snp_1<-snp_1[snp_1[,1] %in% snp_loc[,1],]
expr_1<-expr_1[expr_1[,1] %in% expr_loc[,1],]

colnames(snp_1)[1]<-"id"
colnames(expr_1)[1]<-"id"


snp_random=paste(snpInfo,".SampleSize-",sampleSize,".",prefix,".",chr,".txt",sep="")
expr_random=paste(exprInfo,".SampleSize-",sampleSize,".",prefix,".",chr,".txt",sep="")


write.table(snp_1,snp_random,sep="\t",quote =F, row.names=F,col.names =T)
write.table(expr_1,expr_random,sep="\t",quote =F, row.names=F,col.names =T)


qqplotFile=paste("qqplot.matrixEQtlCis.SampleSize-",sampleSize,".",prefix,".",chr,".pdf",sep="")
outputFileCis=paste("cis-qtl.matrixEQtlAll.SampleSize-",sampleSize,".",prefix,".",chr,".txt",sep="")
outputFileTrans=paste("trans-qtl.matrixEQtlAll.SampleSize-",sampleSize,".",prefix,".",chr,".txt",sep="")
meQTL.random<-matrixEQtlCis(snpInfo=snp_random, exprInfo=expr_random, covarInfo=cov_random, qqplotFile=qqplotFile, snpLoc=snp_chr, exprLoc=expr_chr, outputFileCis=outputFileCis, outputFileTrans=outputFileTrans,pvalueCis=1, pvalueTrans=1e-80, pvalue.hist=FALSE)

file.remove(snp_random)
file.remove(expr_random)
file.remove(snp_chr)
file.remove(expr_chr)
if(!is.null(covarInfo)){
	file.remove(cov_random)
}

