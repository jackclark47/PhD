#Script to concactenate every fasta file for which we have igr data into a single file

##Plan
#Get list of igr loci using the spreadsheet of igr data
#move into the directory containing the entire neisseria pubmlst locus database
#extract the first allele of each locus contained in the list of igr data
#append each locus into a list of sequences
#write that list to a file

setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/")

#Read igr data to get list of loci
loci <- read.xlsx("igr_new_29112023.xlsx")[,1]

#move to directory containing sequence data
setwd("~/Documents/PhD/PhD/N_meningitidis_loci")

#build list of sequences
sequences <- list()
removables <- c()
for(i in 1:length(loci)){
  filename <- paste(loci[i], '.fas', sep = '')
  if(length(list.files()[which(list.files() == filename)]) >0){
    seq <- read.fasta(list.files()[which(list.files() == filename)], as.string = TRUE)[[1]][1]
    sequences <- append(sequences, seq)
  }
  else{
    removables <- c(removables, i)
  }
}

loci <- loci[-c(removables)]

#write to file
write.fasta(sequences, names = loci, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/all_coding_first_alleles.fasta')
             