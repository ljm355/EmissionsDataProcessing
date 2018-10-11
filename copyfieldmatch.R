#run in command line as follows:
#Rscript b:/HestiaSpatialAllocation.R "B:/baltimore/shared.csv" "B:/baltimore/output.csv"

library("foreign")
rm(list=ls())
sourceFile = ""
destFile = ""
srcID = ""
destID = ""
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 1) {
  
  sourceFile = args[1]
  destFile = args[2]
  srcID = args[3]
  destID = args[4]
}
src_tb <- read.dbf(sourceFile)
dest_tb <- read.dbf(destFile)

if(srcID == "FID")
{
  src_tb[,srcID] <- 1:nrow(src_tb)
}
src_tb[,srcID] <- src_tb[,srcID] - 1

matchtb <- match(dest_tb[,destID],src_tb[,srcID])
for(n in 5:length(args)) 
{
  if(args[n] %in% colnames(src_tb)){
    
     dest_tb[,args[n]]<- src_tb[matchtb,args[n]]

  }
}
write.dbf(dest_tb, destFile)