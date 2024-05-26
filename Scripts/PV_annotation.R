#Script to automate identification of PV genes based on phasomeit data. 

#Approach
#Load in xlsx file of phasomeit summary, and .gbks for all 10 isolates
#For each locus in the summary, identify its name and search the corresponding isolate's .gbk
#for that locus.
#Extract the coordinates of the locus from the .gbk
#Extract the sequence using those coordinates from the .fasta file of the isolate
#Blast the sequence against my pubmlst database to find the closest matching locus with a pubmlst id
#assign that id to the sequence
#save the sequence to a fasta file
#check all sequences for that locus to see they all have the same best match
#align those sequences. 
#check the repeat, its length, and on/off state
#extract igrs if repeat is predicted in igr.
library(openxlsx)
library(gggenomes)
library(rBLAST)
library(tidyselect)
test <- read_gbk("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/Annotations/27509/27509.gbff")
#Load data
PV <- read.xlsx("~/Documents/PhD/PhD/PhasomeIt_data/phasomeit_out.xlsx", sheet = 4)[-c(2),]
#load gbk file when looping for that isolate
isolate <- '27509'
locus <- 'JJDFEL_00005'

#function that takes a locus name and isolate as input and extracts the relevant dna sequence
seq_extractor <- function(isolate, locus){
  sequence <- read.fasta(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/Annotations/", isolate, '/', isolate, '.ffn', sep = ''), as.string = T)
  sequence <- sequence[which(substr(names(sequence), nchar(names(sequence))-4,nchar(names(sequence))) == locus)]
  return(as.character(sequence))
}

#Function that takes a sequence and blasts it against my blast db to find the top hit and assign its id
identifier <- function(sequence, blast_db){
  sequence <- DNAStringSet(sequence)
  sequence <- Biostrings::translate(sequence)
  if(length(sequence) == 0){
    print('seq length is 0')
    return(FALSE)
  }
  res <- predict(blast_db, sequence, type = 'blastp')[1,]
  if(is.na(res$evalue)){
    print('no e value')
    return(FALSE)
  }
  if(res$evalue < 5e-10){
    id <- res$sseqid
    
    return(id)
  }
  print('e value too high')
  print(res$evalue)
  print(res)
  return(FALSE)
}
res
seq <- seq_extractor('27509', '04170')
seq
identifier(seq, blast_db)

seq <- seq_extractor('28287', '05855')
seq
identifier(seq, blast_db)


main <- function(PV_data){
  blast_db <- blast(db="~/Documents/PhD/PhD/N_meningitidis_loci/translations/BLAST_db/coding.fsa", type = 'blastp')
  
  for(i in 1:nrow(PV)){
    ids <- c()
    for(j in 10:17){
      isolate <- substr(colnames(PV)[j], start =1, stop = nchar(colnames(PV)[j])-4)
      locus <- PV[i,j]
      locus <- substr(locus, 6, nchar(locus))
      if(is.na(locus)){
        next
      }
      entry <- PV[i,2]
      print(locus)
      sequence <- seq_extractor(isolate, locus)
      if(identifier(sequence, blast_db) == FALSE | length(identifier(sequence, blast_db)) == 0){
        next
      }
      name <- identifier(sequence, blast_db)
      ids <- c(ids, name)
      write.fasta(sequences = sequence, names = name, file.out = paste('~/Documents/PhD/PhD/PhasomeIt_data/PV_loci_fastas/', entry, '_', isolate, '_', locus, '.fasta', sep = ''))
      
    }
    ids <- unique(ids)
    if(length(ids) > 0){
      print(ids)
      print(length(ids))
      PV$PubMLST_id[i] <- ids 
    }
  }
  return(PV)
}

test <- main(PV)

#Now write this sheet to the PV xlsx
write.xlsx(x = test, file = '~/Documents/PhD/PhD/PhasomeIt_data/PV_IDs2.xlsx')
l <- c('a', 'b')
x <- unlist(l)



