#IGR Isolator. Isolates igrs based on sequences extracted from PubMLST. Separate script elsewhere for extraciting igrs identified through BLAST

#Load libraries
library(seqinr)
library(tidyverse)
library(openxlsx)

setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data")
#Read file of loci with IGRs
IGR <- read.fasta("MenYwithflank.xmfa", as.string = TRUE, whole.header = TRUE)
#Read file of loci without IGRs
Coding <- read.fasta("MenYnoflank.xmfa", as.string = TRUE, whole.header = TRUE)
#Extract IGR sequences


#Combine each dataframe (one per isolate) into a larger one where each locus has the entries for each isolate

#Remove the = signs at the end of every 8 entry
for(i in 1:length(Coding)){
  Coding[[i]][1] <- str_replace_all(Coding[[i]][1], '=', '')
  IGR[[i]][1] <- str_replace_all(IGR[[i]][1], '=', '')
}

#Loop through combined df to produce a matrix of allelic differences for each isolates IGRs. 

#Initialise lists to store a list of genes with and without the flanking regions, plus a list of the locus names
gene_with_flanking <- list()
gene_without_flanking <- list()
locus_names <- list()

#for locus in the data, add the relevant sequence to the relevant list (either coding and IGR, or just coding). Track locus names as well
for(gene in 1:length(IGR)){
  locus <- names(IGR)[gene]
  gene_with_flanking <- append(gene_with_flanking, IGR[[gene]][1])
  gene_without_flanking <- append(gene_without_flanking, Coding[[gene]][1])
  locus_names <- append(locus_names, names(IGR)[gene])
}

#Initilaise lists containing just upstream or just downstream IGRs for given loci
upstream_IGR <- list()
downstream_IGR <- list()

#For locus in the data, find start and end coords of the coding sequence, then subtract those from the IGRs to add just their sequence to the relevant list
for(i in 1:length(Coding)){
  upstream_IGR <- append(upstream_IGR, str_split(IGR[[i]][1], Coding[[i]][1])[[1]][1])
  downstream_IGR <- append(downstream_IGR, str_split(IGR[[i]][1], Coding[[i]][1])[[1]][2])
  #find position for that locus seq in IGR list that matches the start of the same locus in the noIGR
}

#Now need to identify IGRs that are smaller than 1000 bases. If this occurs then all IGRs of a given locus are to be removed.
count = 0
removables <- c()
for(i in 1:length(upstream_IGR)){
  if(nchar(upstream_IGR[i]) < 500 | nchar(downstream_IGR[i]) < 500){
    count = count + 1
    print(i)
    index <- 1+ floor((i-1)/8)*8
    removables <- c(removables, index, index+1, index+2, index+3, index+4, index+5, index+6, index+7)
  }
}

count #10640 out of 24800 entries are smaller than 1000. 9769 are smaller than 500 in either the down or upstream region
removables
removables <- unique(removables)
length(removables)
length(locus_names)
#Remove loci which have at least one entry included in removables
remainingloci <- locus_names[-c(removables)]
removedloci <- locus_names[c(removables)]
length(remainingloci) #about 50% removed. 
12680/8 #1585 loci need BLASTing
12120/8 #1515 loci are okay. 

###Process removed loci - remove the duplicates (as 8 for each locus) and make sure sequence names are same as in PubMLST
removedloci <- removedloci[seq(1, length(removedloci), 8)]
removedlocivec <- c()
for(i in 1:length(removedloci)){
  removedloci[i] <- gsub(".*\\s", "", removedloci[i])
  removedlocivec <- c(removedlocivec, removedloci[[i]]) #Convert list to vector
}

#Save removed loci as text file, one entry per line
write(removedlocivec, "IGR_via_BLAST_loci.txt")

#Now remove the actual IGRs
upstream_IGR <- upstream_IGR[-c(removables)]
downstream_IGR <- downstream_IGR[-c(removables)]

#Now shorten the remaining IGRs to 500 bases
for(i in 1:length(upstream_IGR)){
  upstream_IGR[[i]] <- substr(upstream_IGR[[i]], start = (nchar(upstream_IGR[[i]])-499), stop = nchar(upstream_IGR[[i]]))
  downstream_IGR[[i]] <- substr(downstream_IGR[[i]], start = 1, stop = 500)
}

#Write the igrs to files
write.fasta(upstream_IGR, names = remainingloci, file.out = '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/PubMLST_upstream_igr.xmfa')
write.fasta(downstream_IGR, names = remainingloci, file.out = '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/PubMLST_downstream_igr.xmfa')

length(remainingloci)
length(locus_names)
length(upstream_IGR)

#Now compare each group of 8 entries to see if the sequences are identical. By doing so, assign allele numbers. Store as a spreadsheet
out <- matrix(nrow = length(upstream_IGR)/8, ncol = 19)
colnames(out) <- c("Locus", "27509_down", "27553_down", "28262_down", "28269_down", "28287_down",
                   "53930_down", "53948_down", '53951_down', "27509_up", "27553_up", "28262_up", "28269_up",
                   "28287_up", "53930_up", '53948_up', "53951_up", "Alleles-Down", "Alleles-Up")

out <- as.data.frame(out)
#Add locus information


#Dataframe structure. Return one for each isolate

#Locus    #Upstreamcoords   #Upstreamseq    #Downstreamcoords   #Downstreamseq    #Isolate
#NEIS0001         1-131     GCAACTTCATACA     2549-2871           TTCGCTAACCG
#NEIS0002
#NEIS0003
  
  
#Add allele information for the upstream IGRs 
row = 0
for(i in seq(from = 1, to = length(upstream_IGR), by = 8)){
  row = row + 1
  #add locus info to the dataframe
  out$Locus[row] = gsub(".*\\s", "", remainingloci[i])
  
  #make vector of the 8 sequences at that locus
  all_sequences <- c(upstream_IGR[[i]], upstream_IGR[[i+1]], upstream_IGR[[i+2]],upstream_IGR[[i+3]],
                 upstream_IGR[[i+4]],upstream_IGR[[i+5]],upstream_IGR[[i+6]],upstream_IGR[[i+7]])
  #Check how many unique sequences there are
  sequences <- unique(all_sequences)
  #Add the number of unique sequences to the Alleles-up column
  out$`Alleles-Up`[row] <- length(sequences)
  #Track the unique sequences for each locus so that they can be mapped back to the individual isolates and numbered individually
  for(j in 1:length(all_sequences)){
    for(k in 1:length(sequences)){
      if(all_sequences[j] == sequences[k]){
        out[row,(j+9)] <- k
      }
    }
  }
  #repeat for downstream
  #make vector of the 8 sequences at that locus
  all_sequences <- c(downstream_IGR[[i]], downstream_IGR[[i+1]], downstream_IGR[[i+2]], downstream_IGR[[i+3]],
                     downstream_IGR[[i+4]], downstream_IGR[[i+5]], downstream_IGR[[i+6]], downstream_IGR[[i+7]])
  #Check how many unique sequences there are
  sequences <- unique(all_sequences)
  #Add the number of unique sequences to the Alleles-up column
  out$`Alleles-Down`[row] <- length(sequences)
  #Track the unique sequences for each locus so that they can be mapped back to the individual isolates and numbered individually
  for(j in 1:length(all_sequences)){
    for(k in 1:length(sequences)){
      if(all_sequences[j] == sequences[k]){
        out[row,(j+1)] <- k
      }
    }
  }
}

out$Locus[0:100]
out2 = out

#Now compare my results to those in the existing igr spreadsheet
old_igr <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/NEIS_IGR_Alleles-cdb-21092022.xlsx", sheet = 4)
old_igr <- old_igr[order(old_igr$X1),]
out$Locus <- as.character(out$Locus)
out <- out[order(out$Locus),]


colnames(old_igr)[1] <- "Locus"

#For loci contained in out, check if they are the same allele distribution as in old_igr
truecount = 0
falsecount = 0
errorloci = 0
for(locus in 1:length(out$Locus)){
  if(out$Locus[locus] %in% old_igr$Locus){
    row <- match(out$Locus[locus], old_igr$Locus)
    runfalse = 0
    for(j in 2:length(out[locus,])){
      if(out[locus,j] == old_igr[row, j]){
        truecount = truecount + 1
        print(row)
        print(j)
        print(out$Locus[locus])
        print(old_igr$Locus[row])
        print(out[locus,j])
        print(old_igr[row, j])
        next
      }
      if(out[locus,j] != old_igr[row, j]){
        runfalse = runfalse +1
        falsecount = falsecount +1
      }
    }
    if(runfalse > 0){
      errorloci = errorloci+1
    }
  }
}
errorloci #Count of LOCI with at least one error
truecount #raw count of CELLS that match in both dataset
falsecount #raw count of CELLS that do not match in both datasets
(truecount+falsecount)/18

#summary so far
#upstream and downstream igrs seem inverted in the original data. explains the earlier results where there was a surprising uptick in the downstream igrs significance compared to the upstream
#Would have thought itd be the other round, and indeed it is

#The 1000 flank approach gets us roughly 50% of the loci that are contained wihtin the original dataset (but 1/3 of the originals are all NAs)
#Roughyl 50% do not have at least 500 bases in both flanks and are removed
#Of those over 500 bp long, nearly all are 1000 bps long.

#So of the remaining 1515 loci, all are contained in the original igr dataset
#but 872, 57% of the 1515, do not match the data of the original. not sure why that is.
#can pull the sequences of some of the mismatches and check if theres some difference there.
#Only 602, or 40% have an error in the total number of alleles so i think this is an issue of the order in which allele numbers are assigned

#this loop checks cols 18 and 19 (the sum totals of up and downstream alleles to see if they match across datasets)
truecount = 0
falsecount = 0
errorloci = 0
loci_to_check2 <- c()
for(locus in 1:length(out$Locus)){
  if(out$Locus[locus] %in% old_igr$Locus){
    row <- match(out$Locus[locus], old_igr$Locus)
    runfalse = 0
    for(j in 18:19){
      if(out[locus,j] == old_igr[row, j]){
        truecount = truecount + 1
        # print(row)
        # print(j)
        # print(out$Locus[locus])
        # print(old_igr$Locus[row])
        # print(out[locus,j])
        # print(old_igr[row, j])
        next
      }
      if(out[locus,j] != old_igr[row, j]){
        runfalse = runfalse +1
        falsecount = falsecount +1
      }
    }
    if(runfalse > 0){
      errorloci = errorloci+1
      loci_to_check2 <- c(loci_to_check2, out$Locus[locus])
    }
  }
}
errorloci #414 loci dont match the number of alleles present.
truecount
falsecount
loci_to_check2
(truecount+falsecount)/2 #equals 1496, the number of loci we have that are present in both datasets (1515-1496 = 19 loci in my dataset not present in sheet 4 of the excel file)

###Also need to check the matches between my upstream and old_igr downstream and vice versa
truecount = 0
falsecount = 0
errorloci = 0
for(locus in 1:length(out$Locus)){
  if(out$Locus[locus] %in% old_igr$Locus){
    row <- match(out$Locus[locus], old_igr$Locus)
    runfalse = 0
    for(j in 10:17){
      if(out[locus,j] == old_igr[row, j]){
        truecount = truecount + 1
        print(row)
        print(j)
        print(out$Locus[locus])
        print(old_igr$Locus[row])
        print(out[locus,j])
        print(old_igr[row, j])
        next
      }
      if(out[locus,j] != old_igr[row, j]){
        runfalse = runfalse +1
        falsecount = falsecount +1
      }
    }
    if(runfalse > 0){
      errorloci = errorloci+1
    }
  }
}
errorloci


#Next is to get a list of erronuous loci
count = 0
for(i in 1:length(loci_to_check)){
  old_index <- match(loci_to_check[i], old_igr$Locus)
  new_index <- match(loci_to_check[i], out$Locus)
  if(isTRUE(out[new_index,] == old_igr[old_index,]) == FALSE){
    count = count + 1 
    #loop through the conditional output and track which column has most errors
    output <- out[new_index,] == old_igr[old_index,]
    for(j in 1:length(output)){
      #Tally columns where the falses appear
    }
  }
  
}
loci_to_check[i]
out[new_index, ] == old_igr[old_index,]
count


x <- loci_to_check2[130]
old_igr[old_igr$Locus==x,]
out[out$Locus==x,]


output <- old_igr[old_igr$Locus==x,] == out[out$Locus==x,]
output

#Tyhis is checking the two columns which seem to inverse each other when compared between my dataset and the original.
tally = 0
for(i in 1:length(loci_to_check2)){
  x <- loci_to_check2[i]
  output <- old_igr[old_igr$Locus==x,] == out[out$Locus==x,]
  count = 0
  for(j in 1:length(output)){
    if(j == 4| j == 12){
      if(output[j] == FALSE){
        count = count + 1
      }
    }
  }
  if(count > 1){
    tally = tally + 1
  }
}
tally
351+44 = 395 #number of loci which have this inversion between two colimns between the two datasets. 
#351 show false in both 
#could run an alignment of some of these falso sequences to see what they look like. 
