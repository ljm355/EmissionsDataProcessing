#run in command line as follows:
#Rscript b:/HestiaSpatialAllocation.R "B:/baltimore/shared.csv" "B:/baltimore/output.csv"

library("foreign")
rm(list=ls())
sourceFile = ""
destFile = ""
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 1) {
  
  sourceFile = args[1]
  destFile = args[2]
}
src_tb <- read.dbf(sourceFile)
dest_tb = data.frame(src_tb[])
created = FALSE
for(n in 3:length(args)) 
{
  if(args[n] %in% colnames(src_tb) && !created){
    dest_tb <- data.frame(src_tb[,args[n]])
    colnames(dest_tb)[1]<- args[n]  
    created = TRUE
  }
  
  else if(args[n] %in% colnames(src_tb) && created){
    dest_tb[,args[n]]<- src_tb[,args[n]]
  }
}
write.dbf(dest_tb, destFile)