#This script annotates phasomeit outputs by blasting
#hits against a local short copy of the pubmlst database
#containing only one allele of each locus

#libraries
library(gggenomes)

MC58 <- read_gbk("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/MC58.gbff")
MC58fas <- read.fasta("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/MC58.fasta", as.string = T)



get_seq <- function(gene){
  start <- MC58$start[which(MC58$feat_id == gene)]
  stop <- MC58$end[which(MC58$feat_id == gene)]
  strand <- MC58$strand[which(MC58$feat_id == gene)]

  seq <- substr(MC58fas, start = start, stop = stop)
  if(strand == '-'){
    seq <- reverseComplement(DNAString(seq))
  }
  
  return(as.character(seq))
}

gene <- 'gene-NMB0442'
get_seq(gene)
