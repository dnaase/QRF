# downsample.R
# Aug 6, 2015
# 11:32:47 AM
# 
# Author: yaping
###############################################################################

covarInfo=NULL

for (e in commandArgs(TRUE)) {
	ta = strsplit(e,"=",fixed=TRUE)
	if(! is.na(ta[[1]][2])) {
		if(ta[[1]][1] == "wd"){
			wd<-ta[[1]][2] 
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
		if(ta[[1]][1] == "covarInfo"){
			covarInfo<-ta[[1]][2]
		}
	}
}

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


snp_1<-snp[,colnames(snp) %in% random_1_names]
expr_1<-expr[,colnames(expr) %in% random_1_names]
if(!is.null(covarInfo)){
	cov_1<-cov[,colnames(cov) %in% random_1_names]
	cov_1<-cbind(cov[,1],cov_1)
	colnames(cov_1)[1]<-"id"
	cov_random=paste(covarInfo,".SampleSize-",sampleSize,".txt",sep="")
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


snp_random=paste(snpInfo,".SampleSize-",sampleSize,".txt",sep="")
expr_random=paste(exprInfo,".SampleSize-",sampleSize,".txt",sep="")


write.table(snp_1,snp_random,sep="\t",quote =F, row.names=F,col.names =T)
write.table(expr_1,expr_random,sep="\t",quote =F, row.names=F,col.names =T)

