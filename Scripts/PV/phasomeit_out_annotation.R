#This script annotates phasomeit outputs by blasting
#hits against a local short copy of the pubmlst database
#containing only one allele of each locus

#libraries
library(gggenomes)

MC58 <- read_gbk("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/MC58.gbff")
MC58fas <- read.fasta("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/MC58.fasta", as.string = T)

N222 <- read_gbk("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/Annotations/N222.2/N222.2.gbff")
N222fas <- read.fasta("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/N222.2.fasta", as.string = T)

get_seq <- function(gene, gbk, fasta){
  start <- gbk$start[which(gbk$feat_id == gene)]
  stop <- gbk$end[which(gbk$feat_id == gene)]
  strand <- gbk$strand[which(gbk$feat_id == gene)]
  
  seq <- substr(fasta, start = start, stop = stop)
  if(strand == '-'){
    seq <- reverseComplement(DNAString(seq))
  }
  
  return(as.character(seq))
}

x <- as.character(2)
 'log2-fold expression change'

gene <- 'gene-NMB1374'
get_seq(gene, MC58, MC58fas)

gene <- 'cds-N222_00206'
get_seq(gene, N222, N222fas)

