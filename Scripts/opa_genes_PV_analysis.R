#This script queries every sequence identified from phasomeit as being a homologue of an opa gene against a blast database
#the blast database is made from the MC58 reference sequences of each of the 4 opa genes
#the top BLAST result is returned and that annotation assigned to that sequence. 

#Load blast database


#Function to extract the sequence for a given iteration
extract_sequence <- function(sequence_name){
  isolate_id <- substr(sequence_name, 1, 2)
  seqnum <- substr(sequence_name, 3, 7)
  conversions <- as.data.frame(matrix(nrow = 8, ncol = 2))
  colnames(conversions) <- c('strings', 'isolates')
  strings_vec <- c('JJ', 'GF', 'LC', 'CA', 'JC', 'LN', 'CG', 'BM')
  isolates_vec <- c('27509', '27553', '28262', '28269', '28287', '53930', '53948', '53951')
  conversions$strings <- strings_vec
  conversions$isolates <- isolates_vec
  
  isolate <- conversions[which(conversions$string %in% isolate_id),2]
  
  #open the .ffn for that isolate
  fastafile <- read.fasta(paste('~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/Annotations/', isolate, '/', isolate, '.ffn', sep = ''), as.string = TRUE)
  #get the sequence with a number in its name matching seqnum
  sequence <- fastafile[which(substr(names(fastafile), 8, 12) %in% seqnum)]
  sequence <- as.character(sequence)
  return(sequence)
}

#Function to run BLAST on translated versions of the sequence - both original and reverse complement, and return the blast results
PV_blast <- function(sequence){
  blast_db <- blast(db="~/Documents/PhD/PhD/PhasomeIt_data/BLAST_db/opa_seqs", type = 'blastn')
  
  sequence <- DNAStringSet(sequence)
  res <- predict(blast_db, sequence, type = 'blastp')[1,]
  return(res)
}

#Function to take the top blast hit, check it has a good Escore, and assign that to the sequence
assign_id <- function(sequence_name, res, dataset){

  #If the percent identity and evalue are low in the top result, discard
  if(is.na(res$pident) | is.na(res$evalue)){
    next
  }
  #Otherwise, assign the sequence the same name as the top blast hit
  if(res$pident > 80 && res$evalue < 0.005){
    dataset[which(dataset$sequence_names %in% sequence_name), 2] <- res$sseqid
  }
  return(dataset)
}


#main
main <- function(dataset){
  for(sequence_name in dataset$sequence_names){
    print(sequence_name)
    sequence <- extract_sequence(sequence_name)
    print('sequence extracted')
    res <- PV_blast(sequence)
    print('blast complete')
    dataset <- assign_id(sequence_name, res, dataset)
    print('id assigned')
  }
  return(dataset)
}

sequence_name_vec <- c('JJ06005', 'JJ03165', 'JJ09760', 'GF04000', 'GF06270', 'GF09760', 'LC04750', 'LC09860',
'CA09760', 'CA09795', 'JC05445', 'JC09760', 'CG09700', 'BM09705')
dataset <- as.data.frame(matrix(nrow=length(sequence_name_vec), ncol =2))
colnames(dataset) <- c('sequence_names', 'annotation')
dataset$sequence_names <- sequence_name_vec
  
d2 <- main(dataset)
