#R script to process BLAST igr hits
#IGR hit fasta files from the BLAST run contain multiple hits. the coding seq for NEIS00001 will have roughly ten hits in its out.fasta file,
#Each of which can be from different isolates
#I want to trim. these files so they only contain the top hit from the coding seqs own isolate.
#so if the coding seq was from isolate 23262, the only sequence in the file should be the top match against the isolate 23262 genome.


setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out")
#Test an individual sample first
data <- read.fasta("tRNA-val_53951.fastaigr.fasta", as.string = TRUE)
dataname <- tools::file_path_sans_ext("tRNA-val_53951.fastaigr.fasta")
#Trim dataname to contain just the isolate id
dataname1 <- gsub('.fastaigr', '', dataname) #Trim end
dataname <- gsub('.*_', '', dataname1) #Trim start


#For each sequence in the file, if its name doesnt contain the isolate id we want, remove it
#Then take the first element in the file, as lower index = higher e value in this file format. 
dataname %in% gsub('_.*', '', names(data)[13])
dataname
names(data)

indices <- c()
for(i in 1:length(data)){
  print(i)
  data_id <- gsub('_.*', '', names(data)[i])
  if(dataname %in% data_id == FALSE){
    next
  }
  else if(dataname %in% data_id){
    indices <- c(indices, i)
    break
  }
}
indices
data <- data[c(indices)][1]
dataname

#write this sequence into a fasta
write.fasta(data, names = names(data), file.out = paste(dataname1, 'filt.fasta', sep = ''))

#now build into a loop
#here, skip files that are 0 bytes in size as these failed blast for some reason. Store a list of these loci and compare later to the loci missing in the igr_RNA data
emptyfiles <- c()
for(i in 1:length(list.files())){
  print(list.files()[i])
  #Skip files 0 bytes in size and save the name in a list
  #Should I remove all other files of that locus? as they wont be comparable anymore
  if(file.size(list.files()[i]) == 0){
    emptyfiles <- c(emptyfiles, list.files()[i])
    next
  }
  
  if(list.files()[i] == 'fasta_filt' | list.files()[i] == 'best_hit_any'){
    next
  }
  
  data <- read.fasta(list.files()[i], as.string = TRUE)
  dataname <- tools::file_path_sans_ext(list.files()[i])
  #Trim dataname to contain just the isolate id
  dataname1 <- gsub('.fastaigr', '', dataname) #Trim end
  dataname <- gsub('.*_', '', dataname1) #Trim start
  
  #for each sequence in the current file, if its name doesnt contain the id of the isolate we want, remove it. Take the first hit
  indices <- c()
  for(i in 1:length(data)){
    data_id <- gsub('_.*', '', names(data)[i])
    if(dataname %in% data_id == FALSE){
      next
    }
    else if(dataname %in% data_id){
      indices <- c(indices, i)
      break
    }
  }
  data <- data[c(indices)][1]
  
  #write this sequence into a fasta
  write.fasta(data, names = names(data), file.out = paste('fasta_filt/', dataname1, 'filt.fasta', sep = ''))
}

length(emptyfiles) + length(list.files('fasta_filt'))
length(list.files())



#Optional bit of code that takes simply the first hit in the blast output and uses that
for(i in 1:length(list.files())){
  
  print(list.files()[i])
  if(file.size(list.files()[i]) == 0){
    next
  }
  
  if(list.files()[i] == 'fasta_filt' | list.files()[i] == 'best_hit_any'){
    next
  }
  
  data <- read.fasta(list.files()[i], as.string = TRUE)
  
  dataname <- tools::file_path_sans_ext(list.files()[i])
  #Trim dataname to contain just the isolate id
  dataname1 <- gsub('.fastaigr', '', dataname) #Trim end
  
  write.fasta(data[[1]][1], names = names(data)[1], file.out = paste('best_hit_any/', dataname1, 'filt.fasta', sep=''))
}

list.files()
x <- read.fasta(list.files()[1], as.string = TRUE)
x[[1]][1]
names(x)[1]


#Now compare loci in the emptyfiles() list with the list of igrs from Chris. 



#Now compare files in fasta_filt/ and see if those correspond to the igrs from Chris




#Make a multifasta .xmfa file of all the seqs in fasta_filt/ 
setwd("fasta_filt")

for( i in )




for(i in 1:length(list.files))

#Processing and allele assignment of these is done in the IGR_extraction script, slightly modified
