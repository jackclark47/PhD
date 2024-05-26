#Script to check the locus fasta files i downloaded are the full list of those missed in the pubMLST igr extraction
#I then implement two methods of picking a reference allele:
  #Simply taking the first allele of each locus
  #Finding the average length of each locus and taking the first allele with a size similar to the mean

#Get the list of loci needed
loci_list <- as.vector(read.csv("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/IGR_via_BLAST_loci.txt", header = FALSE)$V1)

#set wd to the location of the downloaded fasta files
setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/reference_allele_seqs/all_alleles")

#Check that every locus in the directory is contained within the list
removables <- c()
for(i in 1:length(list.files())){
  locus <- substr(list.files()[i],1 , (nchar(list.files()[i])-4))
  if(!(locus %in% loci_list)){
    removables <- c(removables)
  }
}

#Check every locus in the list has a file in the directory
missing <- c()
for(i in 1:length(loci_list)){
  locus <- paste(loci_list[i], ".fas", sep = '')
  if(!(locus %in% list.files())){
    missing <- c(missing, locus)
  }
}
#15 missing, not available on PubMLST to download. Will just ignore those

#Now obtain the first allele of each and write each to the directory above
for(i in 1:length(list.files())){
  allele <- read.fasta(list.files()[i], as.string = TRUE)[1]
  write.fasta(allele, names(allele), file.out=paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/reference_allele_seqs/", list.files()[i], sep = ''), as.string=TRUE)
}

