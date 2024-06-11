#Script to obtain the coordinates of PV genes in N222.1.2 

###Need to subset the PV loci list to contain just those with genic PV tracts


#We first need a list of PV loci
#Loci which are found in MC58 will use that sequence as the query
#Loci not in MC58 will use the sequence from one of the carriage isolates as the query
#Loci in the list which are missing in the directory will be returned as a vector
library(openxlsx)
library(rBLAST)
library(stringr)
library(magrittr)
library(Biostrings)

#Load PV loci list
PV_loci <- unique(read.xlsx('~/Documents/PhD/PhD/PhasomeIt_data/phasomeit_out.xlsx', sheet = 4)[1:79,1])

#Find loci not contained in 'Phasomeit_date/PV_loci/PV_fastas
joined_names <- str_extract(list.files('~/Documents/PhD/PhD/PhasomeIt_data/PV_loci/PV_fastas'), '[:graph:]+(?=.fasta)')
ids <- str_extract(joined_names, '(?<=[:punct:])[:alnum:]+')
ids <- unique(ids)

missing <- c()
for(i in 1:length(PV_loci)){
  if(PV_loci[i] %in% ids){
    next
  }
  missing <- c(missing, PV_loci[i])
}
missing

#Also find loci for which there is a file for the sequence, but the sequence is empty in every isolate
files <- list.files('~/Documents/PhD/PhD/PhasomeIt_data/PV_loci/PV_fastas')
filesizes <- file.size(files)
files <- cbind(files, filesizes)

files <- as.data.frame(files)
temp <- str_extract(files$files, '[:graph:]+(?=.fasta)')
files$files <- str_extract(temp, '(?<=[:punct:])[:alnum:]+')
#Order by id
files <- files[order(files$files),]

missingseq <- c()
for(i in seq(1,length(files$files), by = 9)){
  filesizes <- files$filesizes[i:(i+8)]
  print(filesizes)
  if(all(as.numeric(filesizes) < 30)){
    print(files$files[i])
    missingseq <- c(missingseq, files$files[i])
  }
  print('======')
}
missingseq

filesizes <- c(1,2,3,4,5,6,7,8)
all(filesizes < 7)
x == TRUE
missing
#remove missing loci from PV_loci
indices <- which(PV_loci %in% missing)
PV_loci <- PV_loci[-indices]



#Function that returns TRUE or FALSE depending on whether the PV loci exist in MC58
ref_check <- function(pubmlst_id){
  filepath <- '~/Documents/PhD/PhD/PhasomeIt_data/PV_loci/PV_fastas'
  ref_files <- list.files(filepath)[which(substr(list.files(filepath), 1,3) == '240')]
  ref_ids <- str_extract(ref_files, pattern = '(?<=[:punct:])[:alnum:]+(?=[:punct:])')
  locus <- ref_ids[which(ref_ids == pubmlst_id)]
  index <- which(ref_ids == pubmlst_id)
  if(length(locus) == 1 && file.size(ref_files[index]) > 30){
    return(TRUE)
  }
  return(FALSE)
}

#Function that returns the sequence to be used as BLAST input
get_seq <- function(pubmlst_id){
  #Check if the sequence exists in MC58
  if(ref_check(pubmlst_id)){
    queryseq <- DNAStringSet(unlist(read.fasta(paste('~/Documents/PhD/PhD/PhasomeIt_data/PV_loci/PV_fastas/240_', pubmlst_id, '.fasta', sep = ''), as.string = T)))
    return(queryseq)
  }
  else{
    carriage_files <- list.files('~/Documents/PhD/PhD/PhasomeIt_data/PV_loci/PV_fastas')
    carriage_files <- carriage_files[which(substr(carriage_files, 1,3) != '240')]
    #Find the files in each isolate with the current pubmlst_id
    carriage_ids <- str_extract(carriage_files, pattern = '(?<=[:punct:])[:alnum:]+(?=[:punct:])')
    indices <- which(carriage_ids == pubmlst_id)
    carriage_files <- carriage_files[indices]
    file_sizes <- file.size(carriage_files)
    if(length(file_sizes) == 0){
      warning(paste('No files found for locus', pubmlst_id))
      return(NA)
    }
    
    for(i in 1:length(file_sizes)){
      size <- file_sizes[i]
      if(size > 30){
        file <- carriage_files[i]
        queryseq <- DNAStringSet(unlist(read.fasta(paste('~/Documents/PhD/PhD/PhasomeIt_data/PV_loci/PV_fastas/', file, sep = ''), as.string = T)))
        return(queryseq)
      }
      if(i == length(file_sizes)){
        warning(paste('No sequences found in MC58 or carriage isolates for locus', pubmlst_id))
        return(NA)
      }
    }
  }
}

get_seq('NEIS0400')



#Function that runs BLAST and processes the output
PV_blast <- function(query_seq, blast_db){
  res <- predict(blast_db, query_seq)
  print(res)
    if(nrow(res) == 0){
      print('locus not found in N222.1.2.fasta')
      return(FALSE)
    }
  #If the E value and alignment length arent good enough, return FALSE
  if(res$pident[1] <= 80 | res$evalue[1] >=0.005){
    warning('No good BLAST hit identified')
    return(FALSE)
  }
  
  start <- res$sstart[1]
  stop <- res$send[1]
  return(c(start, stop))
}

main <- function(PV_loci){
  out <- read.xlsx('~/Documents/PhD/PhD/PhasomeIt_data/phasomeit_out.xlsx', sheet = 4)
  blast_db <- blast(db = '~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/N222_database/N222.2.fasta')
  for(i in 1:length(PV_loci)){
    print(PV_loci[i])
    #Obtain a sequence to use a blast query, prioritising MC58 sequences if available
    queryseq <- get_seq(PV_loci[i])
    
    if(is.na(queryseq)){
      print('locus not found')
      next
    }
    
    #Run BLAST with the sequence against N222.1.2 and extract the coordinates of the top hit
    coords <- PV_blast(queryseq, blast_db)
    try(
      if(coords == FALSE){
      next
    },
    silent = TRUE)
    out$N222.1.2_start[which(out$PubMLST_id == PV_loci[i])] <- coords[1]
    out$N222.1.2_stop[which(out$PubMLST_id == PV_loci[i])] <- coords[2]
  }
  return(out)
}

out <- main(PV_loci)
write.xlsx(out, file = '~/Documents/PhD/PhD/PhasomeIt_data/test_coords.xlsx')
