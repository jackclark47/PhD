#Make a feature file of igrs based on N222.1.2 hybrid sequence ready to visualise on the N222.1.2 genome
#type #start  #stop #var  



#two ways to label:
#label each locus itself with up, down, high or low.
#or label each igr depending on if it is up, down, high or low. for a high or low locus, both igrs would be labelled high or low
  #for a down or up locus, only the down or upstream igr would be labelled

library(tidyverse)
library(rBLAST)
#Load igr data
upvar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_all.xlsx', sheet = 7)
downvar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_all.xlsx', sheet = 4) 
lowvar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_all.xlsx', sheet = 6)
highvar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_all.xlsx', sheet = 5)
novar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_all.xlsx', sheet = 2)

#load N222 genbank file
N222_anno <- read_gbk('~/Documents/PhD/RNA_IGR/Isolate_fastas/N222.2.gbk')
N222_seq <- read.fasta('~/Documents/PhD/RNA_IGR/Isolate_fastas/N222.2.fasta',as.string = T)

out <- as.data.frame(matrix(data = NA, nrow = nrow(N222_anno), ncol = 5))
colnames(out) <- c('Name', 'Type', 'Start', 'Stop', 'Strand')
#name is the pubmlst id
#type is the vargroup 
#start is start coord
#stop is stop coord


#Annotate N222 with BLAST search against mini pubmlst database
db <- blast(db="~/Documents/PhD/PhD/N_meningitidis_loci/translations/BLAST_db/coding.fsa", type = 'blastp')

for(i in 1:nrow(N222_anno)){
  
  print(i)
  
  feature <- N222_anno$feat_id[i]
  direction <- N222_anno$strand[i]
  
  if(N222_anno$type[i] != 'CDS'){
    next
  }
  
  start <- N222_anno$start[i]
  end <- N222_anno$end[i]
  
  seq <- DNAStringSet(substr(N222_seq, start, end))
  if(direction == '-'){
    seq <- reverseComplement(seq)
  }
  aa_seq <- translate(seq)
  
  res <- predict(object = db, newdata = aa_seq)[1,]
  pubmlst_id <- res$sseqid
  
  out$Name[i] <- pubmlst_id
  out$Start[i] <- start
  out$Stop[i] <- end
  out$Strand[i] <- direction
  
  if(pubmlst_id %in% upvar$Locus){
    out$Type[i] <- 'UpVar'
  }
  if(pubmlst_id %in% downvar$Locus){
    out$Type[i] <- 'DownVar'
  }
  if(pubmlst_id %in% lowvar$Locus){
    out$Type[i] <- 'LowVar'
  }
  if(pubmlst_id %in% highvar$Locus){
    out$Type[i] <- 'HighVar'
  }
  if(pubmlst_id %in% novar$Locus){
    out$Type[i] <- 'NoVar'
  }
}

out2 <- na.omit(out)
write.csv(out2, '~/Documents/PhD/Asides/igr_map.csv', row.names = F)
#write as genbank file
#write_gff3()

out <- as.data.frame(matrix(data = NA, nrow = nrow(N222_anno), ncol = 5))
colnames(out) <- c('Name', 'Type', 'Start', 'Stop', 'Strand')

#make igr map of just var vs no var
db <- blast(db="~/Documents/PhD/PhD/N_meningitidis_loci/translations/BLAST_db/coding.fsa", type = 'blastp')

for(i in 1:nrow(N222_anno)){
  
  print(i)
  
  feature <- N222_anno$feat_id[i]
  direction <- N222_anno$strand[i]
  
  if(N222_anno$type[i] != 'CDS'){
    next
  }
  
  start <- N222_anno$start[i]
  end <- N222_anno$end[i]
  
  seq <- DNAStringSet(substr(N222_seq, start, end))
  if(direction == '-'){
    seq <- reverseComplement(seq)
  }
  aa_seq <- translate(seq)
  
  res <- predict(object = db, newdata = aa_seq)[1,]
  pubmlst_id <- res$sseqid
  
  out$Name[i] <- pubmlst_id
  out$Start[i] <- start
  out$Stop[i] <- end
  out$Strand[i] <- direction
  
  if(pubmlst_id %in% upvar$Locus){
    out$Type[i] <- 'Var'
  }
  if(pubmlst_id %in% downvar$Locus){
    out$Type[i] <- 'Var'
  }
  if(pubmlst_id %in% lowvar$Locus){
    out$Type[i] <- 'Var'
  }
  if(pubmlst_id %in% highvar$Locus){
    out$Type[i] <- 'Var'
  }
  if(pubmlst_id %in% novar$Locus){
    out$Type[i] <- 'NoVar'
  }
}
out2 <- na.omit(out)
write.csv(out2, '~/Documents/PhD/Asides/igr_map_var.csv', row.names = F)
#for each locus in N222 check which vargroup it belongs to and write that to a feature file