#Script to extract igr sequences from the BLAST output, leaving the coding sequences behind

#Load libraries
library(seqinr)
library(tidyverse)

#First need to remove the '=' symbol at the end of every 8th sequence in the original blast input data.
setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/")
for(i in 1:length(list.files())){
  if(grepl('.fasta', list.files()[i]) == TRUE){
    currentseq <- read.fasta(list.files()[i], as.string = TRUE, whole.header = TRUE)
    currentseq <- str_replace_all(currentseq, '=', '')
    write.fasta(currentseq, names = names(currentseq), file.out = list.files()[i])
  }
}

setwd("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/best_hit_any/")
#Remove any loci that aren't present in all 8 isolates, as well as sequences shorter than 1000 bases in length(500 bases either side)
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
#remove each file with a file name contained in filenames(), if so, move to dir = under_1000
###An error in this code means that it removes files each iteration, which means it skips some files (as when i is removed, what was i+1 is now i. so the next irteration moves to what was i+2 instead)
#better solution is copy all files over then remove form original directory based on name matches between the two, using the destination folder as the object for the loop.
# for(f in 1:length(list.files())){
#   if(list.files()[f] %in% filenames){
#     file.copy(from = paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/fasta_filt/", list.files()[f], sep = ''),
#               to = "~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/fasta_filt/under_1000")
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
# length(filenames) #1132 on run 1
# length(list.files()) #4581 on run 1
# filenames
# #move files with filename in filenames() to dir= not_universal ###this loop is coded the wrong way round. the iterable, f, skips files as theyre removed
# for(f in 1:length(list.files())){
#   if(list.files()[f] %in% filenames){
#     file.copy(from = paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/best_hit_any/", list.files()[f], sep = ''),
#               to = "~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/fasta_filt/not_universal")
#     file.remove(paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/best_hit_any/", list.files()[f], sep = ''))
#   }
# }

length(list.files()) #3449 on run 1. 
4581-1132
3448/8 #431 loci left which are present in all 8 isolates
#down to 2546 sequences to extract igrs from after multiple runs of the above loops
#The above removal loops (under 1000 and non universal) all need to be run half a dozen times before they finally get evertything. not implemented super well. 


#Now to extract IGRs from each sequence.
###Directory of seqs with flanking == "~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences/BLAST_out/fasta_filt"
###Directory of seqs without flanking == "~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_sequences"
#The former were obtained from BLAST results using the latter as input. Now we subtract the latter from the former to find the flanking regions



#Initialise lists to store up and downstream igrs for each locus
upstream_igr <- list()
downstream_igr <- list()
locus_names <- list()

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
write.fasta(upstream_igr, names = locus_names, file.out = '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_upstream_igr_anybesthit_nolengthfilter.xmfa')
write.fasta(downstream_igr, names = locus_names, file.out = '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_downstream_igranybesthit_nolengthfilter.xmfa')

#If restarting the code, can load IGRS from a single xmfa below:
#upstream_igr <- read.fasta('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_upstream_igrtemp.xmfa')
#downstream_igr <- read.fasta('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/BLAST_downstream_igrtemp.xmfa')

#Check every set of 8 igrs are identical
tester <- c()
for(i in 1:length(locus_names)){
  adder <- gsub('_.*', '', locus_names[i])
  tester <- c(tester, adder)
}
length(tester)
length(locus_names)
length(unique(tester))
length(locus_names)/8
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

106+24+4+1+2+2+2+30
18+15+15+15+16+17+17+17
# isolates
# isolates <- unlist(isolates)
# isolates
# table(isolates)
# counts_per_locus
# length(downstream_igr)/8
# length(counts_per_locus)
# table(counts_per_locus)

#remove loci in removables and see how many are left
removables <- unique(removables)
length(locus_names) - length(removables) #16 entries i.e. 2 loci left and they're tRNA genes


unique(gsub('_.*', '', locus_names)) #get names of the remaining loci
#two tRNA loci ;-;


#Create master objects which are unchanged
up_igr_master <- upstream_igr
down_igr_master <- downstream_igr
locus_names_master <- locus_names

#Any locus with all igrs at 500 is to be saved.
#Remove IGRs which are smaller than 500 bases in either direction and all others of that locus, as there will no longer be 8 samples.
removables <- c()
for(i in 1:length(upstream_igr)){
  print(i)
  if(is.na(upstream_igr[i]) | is.na(downstream_igr[i])){
    print('is na')
    index <- 1+ floor((i-1)/8)*8
    removables <- c(removables, index, index+1, index+2, index+3, index+4, index+5, index+6, index+7)
    next
  }
  if(nchar(upstream_igr[i]) < 500 | nchar(downstream_igr[i]) < 500){
    print('too short')
    count = count + 1
    index <- 1+ floor((i-1)/8)*8
    removables <- c(removables, index, index+1, index+2, index+3, index+4, index+5, index+6, index+7)
  }
}

length(unique(removables))
removables
length(upstream_igr)-length(unique(removables))
1504/8 #= 188, same as in flank_size(500)
upstream_igr <- up_igr_master[-removables]
downstream_igr <- down_igr_master[-removables]
locus_names <- locus_names_master[-removables]

#save these to a file
write.fasta(upstream_igr, names = locus_names, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/up_igr_BLAST_500.xmfa')
write.fasta(downstream_igr, names = locus_names, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/down_igr_BLAST_500.xmfa')


#any locus with all igrs between 500-100 is to be saved
upstream_igr <- up_igr_master[unique(removables)]
downstream_igr <- down_igr_master[unique(removables)]
locus_names <- locus_names_master[unique(removables)]

#truncate the sequences here to 100 nts
for(i in 1:length(upstream_igr)){
  upstream_igr[[i]][1] <- substr(upstream_igr[[i]][1], start = (nchar(upstream_igr[[i]][1])-99), stop = nchar(upstream_igr[[i]][1]))
  downstream_igr[[i]][1] <- substr(downstream_igr[[i]][1], start = 1, stop = 100)
}
nchar(up_igr_master[[1]][1])-100
nchar(substr(up_igr_master[[1]][1], 212,311))

###As shown below, most of the shorter sequences are in the downstream igrs. so will keep the upstream ones separate i think rather than exlude those too
nchar(upstream_igr)[1:999]
nchar(upstream_igr)[1000:1944]
nchar(downstream_igr)[1:999]
nchar(downstream_igr)[1000:1944]

removables <- c()
for(i in 1:length(upstream_igr)){
  print(i)
  if(is.na(upstream_igr[[i]][1]) | is.na(downstream_igr[[i]][1])){
    index <- 1+ floor((i-1)/8)*8
    removables <- c(removables, index, index+1, index+2, index+3, index+4, index+5, index+6, index+7)
  }
  if(!is.na(upstream_igr[[i]][1])){
    if(nchar(upstream_igr[[i]][1]) < 100){
      upstream_igr[[i]][1] <- NA
      index <- 1+ floor((i-1)/8)*8
      removables <- c(removables, index, index+1, index+2, index+3, index+4, index+5, index+6, index+7)
    }
  }
  if(!is.na(downstream_igr[[i]][1])){
    if(nchar(downstream_igr[[i]][1]) < 100){
      downstream_igr[[i]][1] <- NA
      index <- 1+ floor((i-1)/8)*8
      removables <- c(removables, index, index+1, index+2, index+3, index+4, index+5, index+6, index+7)
    }
  }
}

#any locus with igrs between 500-100 and one or two missing is to be saved, with missing entries given NA values for allele assignment
write.fasta(upstream_igr, names = locus_names, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/up_igr_BLAST_100.xmfa')
write.fasta(downstream_igr, names = locus_names, '~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/down_igr_BLAST_100.xmfa')


