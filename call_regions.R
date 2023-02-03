#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are 4 arguments: if not, return an error
if (length(args)!=4) {
   stop("Four arguments must be supplied (input and output file names, column indexes of condition 1, column indexes of condition 2\ne.g input.txt output.txt 4,5,6 7,8,9", call.=FALSE)
}

source("SwitchRT.R")

EL<-read.table(args[1],header=TRUE)
set1=as.numeric(strsplit(args[3],",")[[1]])
set2=as.numeric(strsplit(args[4],",")[[1]])

DRs = getDistances(EL, set1, set2)
DRs = getPvalues(DRs)
DRs = getQvalues(DRs)
DRs<-data.frame(na.omit(EL[,c(1,2,3,set1,set2)]),as.data.frame(DRs))
write.table(DRs,args[2],quote=F,sep="\t",row.names=F)

