#Script to translate each reference sequence for every locus to protein and save in /translations

library(Biostrings)
library(seqinr)
#setwd
setwd("~/Documents/PhD/PhD/N_meningitidis_loci/")

seqs <- list.files()
for(i in 1:length(seqs)){
  print(i)
  if(substr(seqs[i], (nchar(seqs[i])-3), nchar(seqs[i]) ) != '.fas'){
    next
  }
  seq <- read.fasta(seqs[i])
  proseq <- seqinr::translate(seq[[1]][1:length(seq[[1]])])
  
  write.fasta(proseq, names = substr(seqs[i], 1, (nchar(seqs[i])-4) ), file.out= paste("~/Documents/PhD/PhD/N_meningitidis_loci/translations/", seqs[i], sep=''))
}




#now make multifasta of all of the protein sequences
setwd("~/Documents/PhD/PhD/N_meningitidis_loci/translations/")

seq <- list()
for(i in 1:length(list.files())){
  seq <- append(seq, read.fasta(list.files()[i], as.string = TRUE))
}

write.fasta(seq, names = names(seq), file.out = "coding.fsa")


