library("foreign")
library("stringr")

dir = "E:/Vulcan/gridPrep_SHP_master/emissions/"
field = "ca11"
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 1) {
  
  dir = args[1]
  field = args[2]
}

files <- list.files(dir, pattern = ".dbf")

for (i in seq_along(files)) {
  
  name = files[i]
  len = nchar(name)
  ext = substr(name, len-3, len)
  if(ext == '.dbf')
  {
    filename = paste(dir, files[i], sep="")
    tb <- read.dbf(filename)
    tb[,field][is.na(tb[,field])]<- 0
    write.csv(sum(tb[,field]), paste(filename, ".csv", sep=""), row.names=FALSE) #for testing only
  }

}

