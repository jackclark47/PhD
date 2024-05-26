###Script to separate multifasta files into individual sequence files
setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data")
input <- read.fasta("IGR_for_BLAST_seqs.fa", as.string=TRUE, whole.header=TRUE)

#First thing is to separate sequences and locus details
sequences <- c()
loci <- c()
for(i in 1:length(input)){
  sequences <- c(sequences, input[[i]][1])
  loci <- c(loci, names(input[i]))
}

#Add a '>' at the start of each locus name
for(i in 1:length(loci)){
  loci[i] <- paste('>', loci[i], sep = '')
}

#Remove the = signs at the end of every 8 entry
for(i in 1:length(Coding)){
  Coding[[i]][1] <- str_replace_all(Coding[[i]][1], '=', '')
  IGR[[i]][1] <- str_replace_all(IGR[[i]][1], '=', '')
}

#Write each sequence locus pair into its own file
setwd("BLAST_sequences")
for(i in 1:length(loci)){

  name_start <- gsub(".*\\s", "", loci[i])
  name_end <- gsub(":.*", "", loci[i])
  name_end <- substr(name_end, 2, nchar(name_end))
  name <- paste(name_start, name_end, sep = '_')
  filename <- paste(name, '.fasta', sep = '')
  file <- c(loci[i], sequences[i])
  #
  write(file, filename)
}

substr(gsub(":.*", "", loci[i]), 2, nchar(gsub(":.", "", loci[i])))
