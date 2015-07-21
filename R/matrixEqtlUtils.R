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


# matrixEqtlUtils.R
# Nov 26, 2014
# 10:18:16 AM
# 
# Author: yaping
###############################################################################

##Sample script: Illumina 450K array probes to melation probes position (hg18)
methy450k<-read.table("GPL13534_Illumina450K_Manifest.ID_position.txt",sep="\t",header=T,stringsAsFactors=F)
methy450k.hg18<-cbind(methy450k[,1],paste("chr",methy450k[,5],sep=""),as.numeric(methy450k[,6])-1,methy450k[,6])
colnames(methy450k.hg18)<-c("geneid","chr","s1","s2")
write.table(methy450k.hg18,"GPL13534_Illumina450K_Manifest.ID_position.hg18.txt",sep="\t",quote =F, row.names=F,col.names =T)

########convert Illumina 450K array methylation file into matrix eQTL format
methyInfo<-read.table("GSE39672_matrix_processed.txt",sep="\t",na.strings = "", header=T)
methy<-as.matrix(methyInfo[,1])
rownames(methy)<-methyInfo[,1]
methy_name<-"id"
for(i in seq(2,length(methyInfo[1,]),by=2)){
	methy<-cbind(methy,methyInfo[,i])
	methy_name<-c(methy_name,colnames(methyInfo)[i])
}
colnames(methy)<-methy_name

#######convert SNP hapmap file into matrix eQTL format

snpInfo.CEU<-read.table("genotypes_chr1_CEU_r27_nr.b36_fwd.txt",sep=" ", comment.char = "~", header=T, stringsAsFactors=F)
snpInfo.YRI<-read.table("genotypes_chr1_YRI_r27_nr.b36_fwd.txt",sep=" ", comment.char = "~", header=T, stringsAsFactors=F)
rownames(snpInfo.CEU)<-snpInfo.CEU[,1]
rownames(snpInfo.YRI)<-snpInfo.YRI[,1]
common<-intersect(rownames(snpInfo.CEU),rownames(snpInfo.YRI))

snpInfo<-cbind(snpInfo.CEU[common,],snpInfo.YRI[common,12:length(snpInfo.YRI[1,])])

snpInfoBinary<-cbind(snpInfo[,1:11],t(apply(snpInfo[,12:length(snpInfo[1,])],1, char2binary)))
#SNP location
snpLoc<-snpInfoBinary[,c(1,3:4)]
colnames(snpLoc)<-c("snp","chr","pos")

snpInfoBinary<-snpInfoBinary[,c(1,12:length(snpInfo[1,]))]
colnames(snpInfoBinary)[1]<-"id"

common<-intersect(colnames(methy),colnames(snpInfoBinary))

methy.matrixqtl<-methy[,common]
snp.matrixqtl<-snpInfoBinary[,common]

write.table(snpLoc,"Hapmap_genotypes_chr1_CEU_YRI_r27_nr.SNP_loc.hg18.txt",sep="\t",quote =F, row.names=F,col.names =T)
write.table(snp.matrixqtl,"Hapmap_genotypes_chr1_CEU_YRI_r27_nr.SNP_info.hg18.txt",sep="\t",quote =F, row.names=F,col.names =T)
write.table(methy.matrixqtl,"GSE39672_matrix_processed.Zhang2014_LCL.hg18.txt",sep="\t",quote =F, row.names=F,col.names =T)


