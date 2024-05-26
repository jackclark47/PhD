#Script to assign alleles to all igrs, from both blast and pubmlst
#The script then creates an excel spreadsheet with this information on
#Then splits into groups of high, low, up and down variation groups. 

#libraries

#setwd
setwd('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs')

#load igrs
pubmlst_up <- read.fasta('PubMLST_upstream_igr.xmfa', whole.header = TRUE, as.string = TRUE)
blast500_up <- read.fasta('up_igr_BLAST_500.xmfa', whole.header = TRUE, as.string = TRUE)
blast100_up <- read.fasta('up_igr_BLAST_100.xmfa', whole.header = TRUE, as.string = TRUE)
pubmlst_down <- read.fasta('PubMLST_downstream_igr.xmfa', whole.header = TRUE, as.string = TRUE)
blast500_down <- read.fasta('down_igr_BLAST_500.xmfa', whole.header = TRUE, as.string = TRUE)
blast100_down <-read.fasta('down_igr_BLAST_100.xmfa', whole.header = TRUE, as.string = TRUE)

####ALLELES DONT MATCH THE OLD DATASET



#change the pubmlst locus names to match the format for the blast names
for(i in 1:length(pubmlst_up)){
  locus <- names(pubmlst_up)[i]
  id <- substr(locus, 1, 5)
  locus <- gsub(pattern = '.*\\s', replacement = '', x = locus)
  names(pubmlst_up)[i] <- paste(locus, '_', id, sep = '')
  locus <- names(pubmlst_down)[i]
  id <- substr(locus, 1, 5)
  locus <- gsub(pattern = '.*\\s', replacement = '', x = locus)
  names(pubmlst_down)[i] <- paste(locus, '_', id, sep = '')
}
names(pubmlst_up)
names(blast100_down)


#28262:291772-291773 + NEIS0025
#NEIS1404_53948 #BLAST
#combine all 3 to get a uniform dataset
igrs_up <- c(pubmlst_up, blast500_up, blast100_up)
igrs_down <- c(pubmlst_down, blast500_down, blast100_down)

length(igrs_up)/8
length(igrs_down)/8
#save to a file
write.fasta(igrs_up, names = names(igrs_up), 'final_igrs_up.xmfa')
write.fasta(igrs_down, names = names(igrs_down), 'final_igrs_down.xmfa')

#extract locus names
locus_names_all <- names(igrs_up)
locus_names <- c()
for(i in 1:length(locus_names_all)){
  #get just the locus name, then use unique() on it
  locus_names <- c(locus_names, substr(locus_names_all[i], 1, (nchar(locus_names_all[i])-6)))
}
length(locus_names)
15568/8
1946-1692 #254 loci with matching names
locus_names <- unique(locus_names)
length(locus_names)
length(igrs_up)/8 #doesnt match number of locus names, theres definiteily some duplicates appearing somewhere. 


#assign alleles for each locus.
out <- matrix(nrow = length(igrs_up)/8, ncol = 21)
colnames(out) <- c("Locus", "Length_down", "Length_up", "27509_down", "27553_down", "28262_down", "28269_down", "28287_down",
                   "53930_down", "53948_down", '53951_down', "27509_up", "27553_up", "28262_up", "28269_up",
                   "28287_up", "53930_up", '53948_up', "53951_up", "Alleles-Down", "Alleles-Up")

out <- as.data.frame(out)
########IMPORTANT - convert empty sequence to NA values!!
#Add allele information for the upstream IGRs 
row = 0
for(i in seq(from = 1, to = length(igrs_up), by = 8)){
  print(i)
  row = row + 1
  #add locus info to the dataframe
  out$Locus[row] = locus_names[row] 
  
  #make vector of the 8 sequences at that locus
  all_sequences <- c(igrs_up[[i]], igrs_up[[i+1]], igrs_up[[i+2]], igrs_up[[i+3]],
                     igrs_up[[i+4]], igrs_up[[i+5]], igrs_up[[i+6]], igrs_up[[i+7]])
  #Check how many unique sequences there are
  sequencesold <- unique(all_sequences)
  #remove nas from sequences
  sequences <- c()
  for(j in 1:length(sequencesold)){
    if(sequencesold[j] != 'na'){
      sequences <- c(sequences, sequencesold[j])
    }
    if(sequencesold[j] == 'na'){
      print(sequencesold[j])
    }
  }
  if(length(sequences) == 0){
    out[row, (12:19)] <- NA
    next
  }
  out$Length_up[row] <- nchar(sequences)[1]
  #Add the number of unique sequences to the Alleles-up column
  out$`Alleles-Up`[row] <- length(sequences)
  #Track the unique sequences for each locus so that they can be mapped back to the individual isolates and numbered individually
  for(j in 1:length(all_sequences)){
    for(k in 1:length(sequences)){ 
      if(all_sequences[j] == sequences[k]){
        out[row,(j+11)] <- k
      }
      if(all_sequences[j] == 'na'){
        out[row,(j+11)] <- NA
      }
    }
  }
  
  
  
  #repeat for downstream
  #make vector of the 8 sequences at that locus
  all_sequences <- c(igrs_down[[i]], igrs_down[[i+1]], igrs_down[[i+2]], igrs_down[[i+3]],
                     igrs_down[[i+4]], igrs_down[[i+5]], igrs_down[[i+6]], igrs_down[[i+7]])
  
  #Check how many unique sequences there are
  sequencesold <- unique(all_sequences)
  #remove nas from sequences
  sequences <- c()
  for(j in 1:length(sequencesold)){
    if(sequencesold[j] != 'na'){
      sequences <- c(sequences, sequencesold[j])
    }
    if(sequencesold[j] == 'na'){
      print(sequencesold[j])
    }
  }
  if(length(sequences) == 0){
    out[row, (4:11)] <- NA
    next
  }
  #Add the number of unique sequences to the Alleles-down column
  out$`Alleles-Down`[row] <- length(sequences)
  out$Length_down[row] <- nchar(sequences)[1]
  #Track the unique sequences for each locus so that they can be mapped back to the individual isolates and numbered individually
  for(j in 1:length(all_sequences)){
    for(k in 1:length(sequences)){
      if(all_sequences[j] == sequences[k]){
        out[row,(j+3)] <- k
      }
      if(all_sequences[j] == 'na'){
        out[row,(j+3)] <- NA
      }
    }
  }
}



#assign alleles for each locus in igrs_down


#find number of different alleles at each locus and write to new column


#combine igrs_up alleles and igrs_down alleles into a single dataframe

#write to a spreadsheet.
pubmlst_down[["NEIS1135_27509"]][1] == pubmlst_down[["NEIS1135_27553"]][1] 
pubmlst_down[["NEIS1135_27509"]][1] == pubmlst_down[["NEIS1135_28262"]][1] 
pubmlst_down[["NEIS1135_27509"]][1] == pubmlst_down[["NEIS1135_28269"]][1] 
pubmlst_down[["NEIS1135_27509"]][1] == pubmlst_down[["NEIS1135_28287"]][1] 
pubmlst_down[["NEIS1135_27509"]][1] == pubmlst_down[["NEIS1135_53930"]][1] 
pubmlst_down[["NEIS1135_27509"]][1] == pubmlst_down[["NEIS1135_53948"]][1] 
pubmlst_down[["NEIS1135_27509"]][1] == pubmlst_down[["NEIS1135_53951"]][1] 

#load old spreadsheet and compare allele counts and variation
old_igr <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/NEIS_IGR_Alleles-cdb-21092022.xlsx", sheet = 4)
table(out$`Alleles-Down`, useNA = 'always')
table(old_igr$`Alleles-Down`, useNA = 'always')
table(out$`Alleles-Up`, useNA = 'always')
table(old_igr$`Alleles-Up`, useNA = 'always')
#1943 isolates in the new data, including 10 NAs
length(old_igr$`Alleles-Down`)
#2108 isolates in the old data, excluding NAs
1165/2108


#Write new dataset to an excel spreadsheet Maybe including a column for the sequence length would be ideal as well
#Colour 
write.xlsx(out, 'igr_new_06102023.xlsx', keepNA = TRUE)
