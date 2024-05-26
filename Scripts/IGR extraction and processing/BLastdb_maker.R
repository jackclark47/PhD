#Script to combine the fasta files of all 8 MenY isolates by including a unique id in the header for all contigs of a given isolate
#needed for creating a local blast database of the 8 isolates. 

setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas")


db_formatter <- function(isolate){
  #Read in data
  data <- read.fasta(paste(isolate, ".fasta", sep = ''), as.string= TRUE, forceDNAtolower = FALSE)
  print(length(data))
  for(i in 1:length(data)){
    current <- names(data)[i]
    new <- paste(isolate, '_', current, sep = '')
    names(data)[i] <- new
  }
  return(data)
}

test <- db_formatter("27509")



out <- c(db_formatter("27509"),
db_formatter("27553"),
db_formatter("28262"),
db_formatter("28269"),
db_formatter("28287"),
db_formatter("53930"),
db_formatter("53948"),
db_formatter("53951"))

write.fasta(out, names=names(out), file.out='MenY_db.fasta', as.string= TRUE)



write.fasta(data, names=names(data), file.out='27509.fasta', as.string=TRUE)
