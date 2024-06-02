#This script converts the .txt output from PubMLST 'export sequences' from a multifasta to one fasta file per sequence.
#Filenames are taken from the header of each fasta sequence

#Load packages
library(seqinr)
library(stringr)
library(magrittr)

main <- function(multifasta, outdir){
  names(multifasta) %<>% truncate_ids()
  for(i in 1:length(multifasta)){
    sequence <- multifasta[i]
    write.fasta(sequences = as.character(multifasta[i]), names = names(multifasta[i]),
                file.out=paste(outdir, '/',  names(sequence), '.fasta', sep=''))
  }
  
}

#Function to extract the isolate id and locus tag from the names of each sequence
truncate_ids <- function(seq_names){
  out <- c()
  pattern = paste('(?<=', str_escape('+'), '[:space:])[:graph:]+', sep = '')
  for(name in seq_names){
    locustag <- str_extract(name, pattern)
    isolate <- substr(name, start = 1, stop = 5)
    newname <- paste(isolate, locustag, sep = '_')
    out <- c(out, newname)
  }
  return(out)
}

main(multifasta = read.fasta('~/Documents/PhD/PhD/PhasomeIt_data/PV_loci/PV_fastas.txt',
                             as.string = TRUE, whole.header = TRUE),
     outdir = '~/Documents/PhD/PhD/PhasomeIt_data/PV_loci/PV_fastas')
