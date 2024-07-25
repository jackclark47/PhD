#This script takes the BLAST outputs (igr+coding) and extracts eight hits from each BLAST.out file, one hit for each isolate. It identifies failed BLAST runs and those missing hits in some isolates
#It then does the same for just the coding hits in the second BLAST search
#Finally it subtracts the coding from the igr+coding sequences for each isolate for each locus to obtain just the igrs, stored seaprate as up and downstream sequences

library(stringr)
library(seqinr)
#setwd to the file containing BLAST outputs
setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/reference_allele_seqs/BLAST_out/")

#PV ids 
PV_genes <- c("NEIS", "NEIS", "NEIS")

#First step is to get a list of loci that had no hits and delete the files for them to make future loops run faster. 
#Keep a list of removed loci
removables <- vector(length = length(table(file.size(list.files()))[1])) #Initialise vector of removables of length equal to the number of files of size 0
count = 0 #init count variable for accessing the correct index of removables
for(i in 1:length(list.files())){
  if(file.size(list.files()[i]) == 0){
    count = count + 1
    removables[count] <- list.files()[i]
  }
}

#write the removed loci to a list
removed_loci <- removables
for(i in 1:length(removed_loci)){
  removed_loci[i] <- substr(removed_loci[i], 1, (nchar(removed_loci[i]) -13 ))
}
  
write.csv(removed_loci, file = 'removed_loci.txt')

#Remove the loci from the directory
for(i in 1:length(removables)){
  target <- which(list.files() == removables[i])
  file.remove(list.files()[target])
}


#Check how many of the remaining files are missing a match against one or more isolates
#Isolate ids

isolates <- c("27509", "27553", "28262", "28269", "28287", "53930", "53948", "53951")
incomplete_loci <- list()
for(i in 1:length(list.files())){
  if(list.files()[i] == "removed_loci.txt"){
    next
  }
  print(i)
  locus <- read.fasta(list.files()[i], as.string = TRUE)
  count = 0
  for(j in 1:8){
    if(!(isolates[j] %in% substr(names(locus), start = 1, stop = 5))){
      count = count + 1
    }
  }
  if(count > 0){
    incomplete_loci <- c(incomplete_loci, list.files()[i])
    names(incomplete_loci)[length(incomplete_loci)] <- count
  }

}
incomplete_loci #28 loci are missing in at least one isolate
#18 are present in only 1 isolate
#2 are present in 6 isolates
#8 are present in 7 isolates
#Will remove the loci present in only 1 isolate manually


#Now from each file the top hit for each isolate needs to be extracted and saved as its own file. 
for(i in 1:length(list.files())){
  if(!('.fasta' %in% substr(list.files()[i], (nchar(list.files()[i]) -5), nchar(list.files()[i])) ) ){
    next
  }
  
  locus <- read.fasta(list.files()[i], as.string = TRUE)
  for(j in 1:8){
    if(!(isolates[j] %in% substr(names(locus), start =1, stop = 5))){
      next
    }
    if(isolates[j] %in% substr(names(locus), start = 1, stop = 5)){
      seq <- locus[which(substr(names(locus), 1, 5) == isolates[j])][1]
      write.fasta(seq, names(seq), paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/igrs_from_BLAST/", substr(list.files()[i], 1, (nchar(list.files()[i]))-13), '_', substr(names(seq),1,5),'.fasta', sep = ''))
    }
    
  }
}


###Repeat for the coding sequences obtained from the BLAST
setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/reference_allele_seqs/BLAST_out_coding/")

#First step is to get a list of loci that had no hits and delete the files for them to make future loops run faster. 
#Keep a list of removed loci
removables <- vector(length = length(table(file.size(list.files()))[1])) #Initialise vector of removables of length equal to the number of files of size 0
count = 0 #init count variable for accessing the correct index of removables
for(i in 1:length(list.files())){
  if(file.size(list.files()[i]) == 0){
    count = count + 1
    removables[count] <- list.files()[i]
  }
}

#write the removed loci to a list
removed_loci <- removables
for(i in 1:length(removed_loci)){
  removed_loci[i] <- substr(removed_loci[i], 1, (nchar(removed_loci[i]) -13 ))
}

write.csv(removed_loci, file = 'removed_loci.txt')

#Remove the loci from the directory
for(i in 1:length(removables)){
  target <- which(list.files() == removables[i])
  file.remove(list.files()[target])
}


#Check how many of the remaining files are missing a match against one or more isolates
#Isolate ids

isolates <- c("27509", "27553", "28262", "28269", "28287", "53930", "53948", "53951")
incomplete_loci <- list()
for(i in 1:length(list.files())){
  if(list.files()[i] == "removed_loci.txt"){
    next
  }
  
  locus <- read.fasta(list.files()[i], as.string = TRUE)
  count = 0
  for(j in 1:8){
    if(!(isolates[j] %in% substr(names(locus), start = 1, stop = 5))){
      count = count + 1
    }
  }
  if(count > 0){
    incomplete_loci <- c(incomplete_loci, list.files()[i])
    names(incomplete_loci)[length(incomplete_loci)] <- count
  }
  
}
incomplete_loci #28 loci are missing in at least one isolate
#18 are present in only 1 isolate
#2 are present in 6 isolates
#8 are present in 7 isolates
#Will remove the loci present in only 1 isolate manually


#Now from each file the top hit for each isolate needs to be extracted and saved as its own file. 
for(i in 1:length(list.files())){
  if(!('.fasta' %in% substr(list.files()[i], (nchar(list.files()[i]) -5), nchar(list.files()[i])) ) ){
    next
  }
  
  locus <- read.fasta(list.files()[i], as.string = TRUE)
  for(j in 1:8){
    if(!(isolates[j] %in% substr(names(locus), start =1, stop = 5))){
      next
    }
    if(isolates[j] %in% substr(names(locus), start = 1, stop = 5)){
      seq <- locus[which(substr(names(locus), 1, 5) == isolates[j])][1]
      write.fasta(seq, names(seq), paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/coding_from_BLAST/", substr(list.files()[i], 1, (nchar(list.files()[i]))-13), '_', substr(names(seq),1,5),'.fasta', sep = ''))
    }
    
  }
}


###Extracting the igrs
#if the list.files() in both dirs match then subtract the coding seq from the igr+coding, saving
#the output in a third file in another directory.
setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/igrs_from_BLAST")
list.files("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/igrs_from_BLAST") == list.files("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/coding_from_BLAST")
length(list.files())
length(list.files("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/coding_from_BLAST"))
list.files()[417]
list.files("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/coding_from_BLAST")[417]


up_igrs <- list()
down_igrs <- list()
for(i in 1:length(list.files())){
  setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/coding_from_BLAST")
  coding <- read.fasta(list.files()[i], as.string = TRUE)
  setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/igrs_from_BLAST")
  flank <- read.fasta(list.files()[i], as.string = TRUE)
  
  up_igr <- str_split(flank, coding[[1]][1])[[1]][1]
  down_igr <- str_split(flank, coding[[1]][1])[[1]][2]
  
  up_igrs <- append(up_igrs, up_igr)
  down_igrs <- append(down_igrs, down_igr)
  
}

#Write the up and down igrs into a multifasta
write.fasta(up_igrs, names = list.files(), file.out = '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_upstream_igrs.xmfa')
write.fasta(down_igrs, names = list.files(), file.out = '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_downstream_igrs.xmfa')

#The loci that are not present in all 8 isolates will be separated into another directory as the following code relies on multiples of 8

loci_list <- list.files("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/reference_allele_seqs/BLAST_out")
for(i in 1:length(loci_list)){
  loci_list[i] <- substr(loci_list[i],1, (nchar(loci_list[i])-13))
}

#Check each locus in loci_list is present 8 times in the igrs_from_blast directory
igr_list <- list.files()
for(i in 1:length(igr_list)){
  igr_list[i] <- substr(igr_list[i],1, (nchar(igr_list[i])-12))
}

igr_list[which(igr_list == 'gyrA')]
loci_list[which(loci_list == 'gyrA')]

for(i in 1:length(loci_list)){
  if(str_count(paste(unlist(igr_list), collapse = ''), loci_list[i]) != 8){
    print(loci_list[i])
  }
}

str_count(paste(unlist(igr_list), collapse = ''), 'NEIS0447') #need to double check manually as some loci are present twice in the form of an igr locus etc so are double counted
#Loci missing 8 copies are:
#c("NEIS0447", "NEIS0602", "NEIS1468", "NEIS1469", "NEIS2370", "NEIS2486", "NEIS2490", "NEIS2562", "NEIS2902", "NEIS3045")
#each of these will be moved to a different directory

#Check every set of 8 igrs are identical
tester <- c()
for(i in 1:length(list.files())){
  adder <- gsub('_.*', '', list.files()[i])
  tester <- c(tester, adder)
}

length(unique(tester))
length(list.files())/8 #TRue if you exclude NG_MAST and NG as these are present in 2 and 3 sets of 8 loci, respectively
nchar(down_igrs) == nchar(up_igrs)


#how many isolates are shorter than 500 nt at each locus?
#counts_per_locus stores the number of isolates per locus that are shorter than 500 nt in the upstream
flank_size <- function(seq_length){
  counts_per_locus <- c()
  isolates <- c()
  for(i in seq(1, length(up_igrs), 8)){
    isolate <- list.files()[i:(i+7)]
    for(j in 1:8){
      isolate[j] <- substr(isolate[j],start = (nchar(isolate[j])-4), stop = nchar(isolate[j]))
    }
    seqlengthsup <- nchar(up_igrs[i:(i+7)])
    seqlengthsdown <- nchar(down_igrs[i:(i+7)])
    count = 0
    for(j in 1:8){
      if(is.na(seqlengthsup[j]) | is.na(seqlengthsdown[j])){
        count = count + 1
        next
      }
      if(seqlengthsup[j] < seq_length | seqlengthsdown[j] < seq_length){
        count = count + 1
        isolates <- c(isolates, isolate[j])
      }
    }
    counts_per_locus <- c(counts_per_locus, count)
  }
  print(table(counts_per_locus))
}

flank_size(500)
flank_size(400)
flank_size(300)
flank_size(200)
flank_size(100)

#Trim all igrs for each locus so that they are all the same length, on an igr by igr basis
count = 0
for(i in 1:length(up_igrs)){
  if(is.na(down_igrs[i])){
    next
  }
  if(nchar(down_igrs[i]) < 100){ #can change size limit or down/up igr to see different things
    count = count + 1
  }
}
count #4592 sequences are 500 nt long in the up_igr, 4281 in the down_igr
#447 up_igr sequences are < 100 nts long, 464 down_igrs

#save up and down igrs in holder objects so i dont have to reload the script if the below code goes wrong
up_igrs2 <- up_igrs
down_igrs2 <- down_igrs

up_igrs <- up_igrs2
down_igrs <- down_igrs2

for(i in seq(1,length(up_igrs), 8)){
  
  #upstream truncation
  octet <- c(up_igrs[i], up_igrs[i+1], up_igrs[i+2], up_igrs[i+3], 
             up_igrs[i+4], up_igrs[i+5], up_igrs[i+6], up_igrs[i+7])
  #order in ascending string length, so first index is the shortest sequence
  octet <- octet[order(nchar(octet))]
  
  #take the length of the smallest isolate for the locus with length greater than 99
  for(j in 1:8){
    if(is.na(octet[j])){
      next
    }
    if(nchar(octet[j]) < 100){
      next
    }
    if(nchar(octet[j]) > 99){
      minlength <- (nchar(octet[j])-1) 
      break
    }
  }
  
  #Now truncate all of the sequences so they are the same length as the shortest. 
  for(j in 1:8){
    if(nchar(up_igrs[[i+j-1]][1]) < 100){ #dont substr those seqs shorter than 100. theyll be removed in a later step. 
      next
    }
    if(nchar(up_igrs[[i+j-1]][1]) > 99){
      up_igrs[[i+j-1]][1] <- substr(up_igrs[[i+j-1]][1], start = (nchar(up_igrs[[i+j-1]][1])-minlength), stop = nchar(up_igrs[[i+j-1]][1]))
    }

  }
  
  #Now repeat for downstream, which is cut from the start to a midpoint, whereas upstream is cut from a midpoint to the end
  
  octet <- c(down_igrs[i], down_igrs[i+1], down_igrs[i+2], down_igrs[i+3], 
             down_igrs[i+4], down_igrs[i+5], down_igrs[i+6], down_igrs[i+7])
  octet <- octet[order(nchar(octet))]
  
  #take the length of the smallest isolate for the locus with length greater than 99
  for(j in 1:8){
    if(is.na(octet[j])){
      next
    }
    if(nchar(octet[j]) < 100){
      next
    }
    if(nchar(octet[j]) > 99){
      minlength <- (nchar(octet[j])-1) 
      break
    }
  }
  
  for(j in 1:8){
    if(nchar(down_igrs[[i+j-1]][1]) < 100){ #dont substr those seqs shorter than 100. theyll be removed in a later step. 
      next
    }
    if(nchar(down_igrs[[i+j-1]][1]) > 99){
      down_igrs[[i+j-1]][1] <- substr(down_igrs[[i+j-1]][1], start = 1, stop = minlength+1) #seqs longer than 99 nts are all truncated to the shortest
    }
  }
  
}

#check locus lengths
nchar(up_igrs)
nchar(down_igrs)
table(nchar(up_igrs))
table(nchar(down_igrs))

#convert any sequence below 100 nts in length to an NA value
for(i in 1:length(up_igrs)){
  if(!is.na(up_igrs[[i]][1])){
    if(nchar(up_igrs[[i]][1]) < 100){
      up_igrs[[i]][1] <- NA
    }
  }
  if(!is.na(down_igrs[[i]][1])){
    if(nchar(down_igrs[[i]][1]) < 100){
      down_igrs[[i]][1] <- NA
    }
  }
}

write.fasta(up_igrs, names = list.files(), '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/up_igr_BLAST.xmfa')
write.fasta(down_igrs, names = list.files(), '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/down_igr_BLAST.xmfa')
length(up_igrs)
length(down_igrs)
length(list.files())

nchar(up_igrs)
nchar(down_igrs)
