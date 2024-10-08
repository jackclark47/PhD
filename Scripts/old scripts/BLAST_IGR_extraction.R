#Script to extract igr sequences from the BLAST output, leaving the coding sequences behind

#Load libraries
library(seqinr)
library(tidyverse)

#First need to remove the '=' symbol at the end of every 8th sequence in the original blast input data.
# setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/")
# for(i in 1:length(list.files())){
#   if(grepl('.fasta', list.files()[i]) == TRUE){
#     currentseq <- read.fasta(list.files()[i], as.string = TRUE, whole.header = TRUE)
#     currentseq <- str_replace_all(currentseq, '=', '')
#     write.fasta(currentseq, names = names(currentseq), file.out = list.files()[i])
#   }
# }

# ###Testing
# list.files()
# substr(flank, start = 1, stop = (nchar(flank)-10))
# substr(seq, start = 1, stop = (nchar(seq)-6))
# 
# list.files()[1]
# flank
# seq
# 
# read.fasta(flank, as.string = TRUE, whole.header = TRUE)[[1]]
# read.fasta(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/", seq, sep = ''), as.string = TRUE, whole.header = TRUE)[[1]]
# 
# flankseq <- read.fasta(flank, as.string = TRUE, whole.header = TRUE)[[1]][1]
# seqseq <- read.fasta(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/", seq, sep = ''), as.string = TRUE, whole.header = TRUE)[[1]][1]
# 
# flankseq
# seqseq
# 
# test <- strsplit(flankseq, seqseq)
# test[[1]][2]
# upstream_igr <- append(upstream_igr, str_split(read.fasta(flank, as.string = TRUE, whole.header = TRUE)[[1]][1], read.fasta(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/", seq, sep = ''), as.string = TRUE, whole.header = TRUE)[[1]][1])[1])
# downstream_igr <- append(downstream_igr, str_split(read.fasta(flank, as.string = TRUE, whole.header = TRUE)[[1]][1], read.fasta(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/", seq, sep = ''), as.string = TRUE, whole.header = TRUE)[[1]][1])[2])

#Now check which isolates failed blast and save to an object
# for(i in ){
#   
# }


setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/fasta_filt/")
#Remove any loci that arent present in all 8 isolates, as well as sequences shorter than 1000 bases in length(500 bases either side)
# filenames <- c()
# for(i in 1:length(list.files())){
#   if(grepl('.fasta', list.files()[i]) == TRUE){
#     print(list.files()[i])
#     sequence <- read.fasta(list.files()[i], as.string = TRUE)[[1]][1]
#     #checking if < 1000 bases
#     if(nchar(sequence) < 1000){
#       print('adding')
#       filenames <- c(filenames, list.files()[i])
#     }
#   }
# }
# length(list.files())
# length(filenames)
# #remove each file with a file name contained in filenames(), if so, move to dir = under_1000
# ###An error in this code means that it removes files each iteration, which means it skips some files (as when i is removed, what was i+1 is now i. so the next irteration moves to what was i+2 instead)
# #better solution is copy all files over then remove form original directory based on name matches between the two, using the destination folder as the object for the loop.
# for(f in 1:length(list.files())){
#   if(list.files()[f] %in% filenames){
#     file.copy(from = paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/fasta_filt/", list.files()[f], sep = ''),
#                 to = "~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/fasta_filt/under_1000")
#     file.remove(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/fasta_filt/", list.files()[f], sep = ''))
#   }
# }
# length(list.files())
# length(filenames)

#checking if not present in all 8 isolates, if so, move to dir = not_universal
#initlaise variables to control skipping loop iterations and storing names of loci that arent present in all 8 osolates
# filenames = c()
# skip_count = 0
# skip = FALSE
# #for each sequence run following
# for(f in 1:(length(list.files())-8)){
#   
#   #if skip is set to true and the skip count is less than 7, skip iterations until skip count reaches 7. increase skipcount by 1
#   if(skip == TRUE && skip_count < 7){
#     skip_count = skip_count + 1
# 
#     next
#   }
#   
#   #if skipcount exceeds 6, set skip to FALSE as we should now be on a new locus. 
#   if(skip_count > 6){
#     skip_count = 0
#     skip = FALSE
#   }
#   
#   #the core of the loop. if the list.files object is a fasta file, combine it and the next 7 files into a vector
#   if(grepl('.fasta', list.files()[f]) == TRUE){
#     octet <- c(list.files()[f], list.files()[f+1], list.files()[f+2], list.files()[f+3], 
#                list.files()[f+4], list.files()[f+5], list.files()[f+6], list.files()[f+7])
# 
#     #remove species ids from each entry in octet, getting just the locus name
#     for(i in 1:8){
#       octet[i] <- gsub('_.*', '', octet[i])
#     }
#     
#     #if the loci arent identical, add the filename of the first entry to filenames, ready for removal later
#     if(length(unique(octet)) != 1){
#       filenames <- c(filenames, list.files()[f])
#     }
#     
#     #if the loci are identical, set skip to TRUE so that the next 7 iterations arent run. these would begin incorporating other loci at the end of octet and give false readings, so need to skip them. also makes it run faster
#     if(length(unique(octet)) == 1){
#       skip = TRUE
#     }
#   }
# }
# 
# length(filenames)
# length(list.files())
# 
# #move files with filename in filenames() to dir= not_universal
# for(f in 1:length(list.files())){
#   if(list.files()[f] %in% filenames){
#     file.copy(from = paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/fasta_filt/", list.files()[f], sep = ''),
#               to = "~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/fasta_filt/not_universal")
#     file.remove(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/fasta_filt/", list.files()[f], sep = ''))
#   }
# }
# 
# length(list.files())
# 4582-1132
# 3448/8
#down to 3448 sequences to extract igrs from after multiple runs of the above loops
#The above removal loops (under 1000 and non universal) all need to be run half a dozen times before they finally get evertything. not implemented super well. 


#Now to extract IGRs from each sequence.
###Directory of seqs with flanking == "~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/fasta_filt"
###Directory of seqs without flanking == "~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences"
#The former were obtained from BLAST results using the latter as input. Now we subtract the latter from the former to find the flanking regions

setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/fasta_filt")

#Initialise lists to store up and downstream igrs for each locus
upstream_igr <- list()
downstream_igr <- list()
locus_names <- list()

nchar(read.fasta(list.files()[1], as.string = TRUE)[[1]][1])
nchar(read.fasta(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/", list.files("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/")[2], sep = ''), as.string = TRUE)[[1]][1])
#for each sequence with flank, run the following:
count = 1
for(i in 1:length(list.files())){
  flank = list.files()[i] #seq with flank
  for(j in count:length(list.files("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences"))){
    seq = list.files("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences")[j] #seq without flank
    print(paste('seq is:', seq))
    print(paste('flank is:', flank))
    print('')
    #if they have the same filename, subtract seq from flank
    if(substr(flank, start = 1, stop = (nchar(flank)-10)) == substr(seq, start = 1, stop = (nchar(seq)-6))){
      upstream_igr <- append(upstream_igr, str_split(read.fasta(flank, as.string = TRUE, whole.header = TRUE)[[1]][1], read.fasta(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/", seq, sep = ''), as.string = TRUE, whole.header = TRUE)[[1]][1])[[1]][1])
      downstream_igr <- append(downstream_igr, str_split(read.fasta(flank, as.string = TRUE, whole.header = TRUE)[[1]][1], read.fasta(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/", seq, sep = ''), as.string = TRUE, whole.header = TRUE)[[1]][1])[[1]][2])
      locus_names <- append(locus_names, as.character(substr(flank, start = 1, stop = (nchar(flank)-10))))
      count = j
      break
    }
  }
}

#save igrs now as the above loop takes an exceptionally long time. Because its looping through millions of permutations
write.fasta(upstream_igr, names = locus_names, file.out = '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_upstream_igr_sameisolate_nolengthfilter.xmfa')
write.fasta(downstream_igr, names = locus_names, file.out = '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_downstream_igr_sameisolate_nolengthfilter.xmfa')

#If restarting the code, can load IGRS from a single xmfa below:
#upstream_igr <- read.fasta('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_upstream_igrtemp.xmfa')
#downstream_igr <- read.fasta('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_downstream_igrtemp.xmfa')

#Check every set of 8 igrs are identical
tester <- c()
for(i in 1:length(locus_names)){
  adder <- gsub('_.*', '', locus_names[i])
  tester <- c(tester, adder)
}

length(unique(tester)) == length(locus_names)/8 #TRUE
nchar(downstream_igr) == nchar(upstream_igr)

#how many isolates are shorter than 500 nt at each locus?
#counts_per_locus stores the number of isolates per locus that are shorter than 500 nt in the upstream
flank_size <- function(seq_length){
  counts_per_locus <- c()
  isolates <- c()
  for(i in seq(1, length(upstream_igr), 8)){
    isolate <- locus_names[i:(i+7)]
    for(j in 1:8){
      isolate[j] <- substr(isolate[j],start = (nchar(isolate[j])-4), stop = nchar(isolate[j]))
    }
    seqlengthsup <- nchar(upstream_igr[i:(i+7)])
    seqlengthsdown <- nchar(downstream_igr[i:(i+7)])
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
  print(table(unlist(isolates)))
}

flank_size(500)
flank_size(400)
flank_size(300)
flank_size(200)
flank_size(100)

unique(gsub('_.*', '', locus_names)) #get names of the remaining loci

#Create master objects which are unchanged
up_igr_master <- upstream_igr
down_igr_master <- downstream_igr
locus_names_master <- locus_names
upstream_igr <- up_igr_master
downstream_igr <- down_igr_master
locus_names <- locus_names_master
#Any locus with all igrs at 500 is to be saved.
#Remove IGRs which are smaller than 500 bases in either direction and all others of that locus, as there will no longer be 8 samples.
removables <- c()
for(i in 1:length(upstream_igr)){
  if(is.na(upstream_igr[i])){
    index <- 1+ floor((i-1)/8)*8
    removables <- c(removables, index, index+1, index+2, index+3, index+4, index+5, index+6, index+7)
    next
  }
  if(nchar(upstream_igr[i]) < 500){
    count = count + 1
    index <- 1+ floor((i-1)/8)*8
    removables <- c(removables, index, index+1, index+2, index+3, index+4, index+5, index+6, index+7)
  }
  #shorten any igrs longer than 500
  if(nchar(upstream_igr[i]) > 500){
    print(i)
    print('truncating up')
    upstream_igr[i] <- substr(upstream_igr[i], start = (nchar(upstream_igr[i])-499), stop = nchar(upstream_igr[i]))
  }
}

for(i in 1:length(downstream_igr)){
  if(is.na(downstream_igr[i])){
    index <- 1+ floor((i-1)/8)*8
    removables <- c(removables, index, index+1, index+2, index+3, index+4, index+5, index+6, index+7)
    next
  }
  
  if(nchar(downstream_igr[i]) < 500){
    count = count + 1
    index <- 1+ floor((i-1)/8)*8
    removables <- c(removables, index, index+1, index+2, index+3, index+4, index+5, index+6, index+7)
  }
  
  if(nchar(downstream_igr[i]) > 500){
    print(i)
    print('trunacting down')
    downstream_igr[i] <- substr(downstream_igr[i], start = 1, stop = 500)
  }
}
nchar(upstream_igr[315])
length(unique(removables))
removables
length(upstream_igr)-length(unique(removables))
16/8 #= 2, same as in flank_size(500)
upstream_igr <- up_igr_master[-removables]
downstream_igr <- down_igr_master[-removables]
locus_names <- locus_names_master[-removables]

#save these to a file
write.fasta(upstream_igr, names = locus_names, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/up_igr_BLAST_500.xmfa')
write.fasta(downstream_igr, names = locus_names, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/down_igr_BLAST_500.xmfa')
removables_master <- removables
removables <- removables_master
#any locus with all igrs between 500-100 is to be saved
upstream_igr <- up_igr_master[unique(removables)]
downstream_igr <- down_igr_master[unique(removables)]
locus_names <- locus_names_master[unique(removables)]


###write code here that looks at each 8 sequences in upstream_igr, finds the shortest of the 8 that is >99 nt, and truncates the other 7 sequences to that length
#seqs < 100 nt are assigned NA
for(i in seq(1, length(upstream_igr), 8)){
  #for upstream
  octet <- c(upstream_igr[i], upstream_igr[i+1], upstream_igr[i+2], upstream_igr[i+3],
             upstream_igr[i+4], upstream_igr[i+5], upstream_igr[i+6], upstream_igr[i+7])
  octet <- octet[order(nchar(octet))]
  for(j in 1:8){
    if(is.na(nchar(octet[j]))){
      next
    }
    if(nchar(octet[j]) > 99){
      minlength <- (nchar(octet[j])-1)
      break
    }
  }
  for(j in 1:8){
    upstream_igr[[i+j-1]][1] <- substr(upstream_igr[[i+j-1]][1], start = (nchar(upstream_igr[[i+j-1]][1])-minlength), stop = nchar(upstream_igr[[i+j-1]][1]))
  }
  #for downstream
  octet <- c(downstream_igr[i], downstream_igr[i+1], downstream_igr[i+2], downstream_igr[i+3],
             downstream_igr[i+4], downstream_igr[i+5], downstream_igr[i+6], downstream_igr[i+7])
  octet <- octet[order(nchar(octet))]
  for(j in 1:8){
    if(is.na(nchar(octet[j]))){
      next
    }
    if(nchar(octet[j]) > 99){
      minlength <- (nchar(octet[j])-1)
      break
    }
  }
  for(j in 1:8){
    downstream_igr[[i+j-1]][1] <- substr(downstream_igr[[i+j-1]][1], start = 1, stop = minlength+1)
  }
}


###As shown below, most of the shorter sequences are in the downstream igrs. so will keep the upstream ones separate i think rather than exlude those too
nchar(upstream_igr)[1:999]
nchar(upstream_igr)[1000:1944]
nchar(downstream_igr)[1:999]
nchar(downstream_igr)[1000:1944]

for(i in 1:length(upstream_igr)){
  print(i)
  if(!is.na(upstream_igr[[i]][1])){
    if(nchar(upstream_igr[[i]][1]) < 100){
      upstream_igr[[i]][1] <- NA
    }
  }
  if(!is.na(downstream_igr[[i]][1])){
    if(nchar(downstream_igr[[i]][1]) < 100){
      downstream_igr[[i]][1] <- NA
    }
  }
}

# #Remove loci with all NA values
# removables <- c()
# for(i in seq(1,length(upstream_igr), 8)){
#   octet <- c(upstream_igr[i], upstream_igr[i+1], upstream_igr[i+2], upstream_igr[i+3],
#              upstream_igr[i+4], upstream_igr[i+5], upstream_igr[i+6], upstream_igr[i+7])
#   if(length(unique(octet)) == 1 && is.na(unique(octet))){
#     removables <- c(removables, i, i+1, i+2, i+3, i+4, i+5, i+6, i+7)
#   }
# }
# upstream_igr <- upstream_igr[-removables]
# locus_names_up <- locus_names[-removables]
# 
# removables <- c()
# for(i in seq(1,length(downstream_igr), 8)){
#   octet <- c(downstream_igr[i], downstream_igr[i+1], downstream_igr[i+2], downstream_igr[i+3],
#              downstream_igr[i+4], downstream_igr[i+5], downstream_igr[i+6], downstream_igr[i+7])
#   if(length(unique(octet)) == 1 && is.na(unique(octet))){
#     removables <- c(removables, i, i+1, i+2, i+3, i+4, i+5, i+6, i+7)
#   }
# }
# downstream_igr <- downstream_igr[-removables]
# locus_names_down <- locus_names[-removables]

#any locus with igrs between 500-100 and one or two missing is to be saved, with missing entries given NA values for allele assignment
write.fasta(upstream_igr, names = locus_names, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/up_igr_BLAST_100.xmfa')
write.fasta(downstream_igr, names = locus_names, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/down_igr_BLAST_100.xmfa')
length(upstream_igr)
length(downstream_igr)
length(locus_names)

