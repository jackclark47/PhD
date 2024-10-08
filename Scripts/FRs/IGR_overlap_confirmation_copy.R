#Script to investigate igr variation in loci appearing as hits in significant_loci_matches.xlsx
#IGRs are stored in the final_igrs_[up][down].xlsx files.
#For a given locus:
#Extract the igr seqs for that locus across all 8 isolates.
#Align those sequences so the variable regions can be observed
#Save the alignment to a file
#BLAST one of the up and one of the down igrs against the hybrid reference (N222.2.gbk) to see if it overlaps with a coding region

#Libraries
library(openxlsx)
library(stringr)
library(seqinr)
library(msa)
library(Biostrings)
library(DECIPHER)
library(rBLAST)
library(gggenomes)

#Files
varhits <- hits[which(hits$Group != "NoVar"),]
#only want CDS features from the gbk

#Function to extract all the igrs of a given locus
extractor <- function(locus){
  pattern <- paste('^', locus, '_', sep = "")
  up_locus <- unlist(up_igrs[str_which(names(up_igrs), pattern = pattern)])
  down_locus <- unlist(down_igrs[str_which(names(down_igrs), pattern = pattern)])
  
  #Convert 'na' to NA
  up_locus[which(up_locus == 'na')] <- NA
  down_locus[which(down_locus == 'na')] <- NA
  
  return(list(up_locus, down_locus))
}
igrs <- extractor("NEIS0320")
igrs
#Function to align the igrs of a locus and write to a file
align <- function(igrs, locus){
  #upstream
  file_name <- paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/significant_matches_alignments/", locus, 'up.align', sep = '')
  alignment <- msa(na.omit(unlist(igrs[1])), type = 'dna')
  write.phylip(alignment, file_name)
  alignment <- DNAStringSet(alignment)
  up_cons <- ConsensusSequence(alignment, minInformation = 0.9)
  alignmentup <- DNAStringSet(as.character(alignment))
    
  #downstream
  file_name <- paste("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/significant_matches_alignments/", locus, 'down.align', sep = '')
  alignment <- msa(na.omit(unlist(igrs[2])), type = 'dna')
  write.phylip(alignment, file_name)
  alignment <- DNAStringSet(alignment)
  down_cons <- ConsensusSequence(alignment, minInformation = 0.9)
  alignment <- DNAStringSet(as.character(alignment))
  #extract consensus sequences and the alignments
  return(list(c(up_cons, alignmentup), c(down_cons, alignment)))
}
alignment = align(igrs, "NEIS0361")
alignment
alignment[3:10]
nchar(igrs[[1]][1])


feat_overlap(1792, 2895, c(1771, 2213), 'forward', c(0,0))
feat_overlap <- function(feature_start, feature_stop, hitcoords, direction, overlaps){
  
  overlap_for <- overlaps[1]
  overlap_rev <- overlaps[2]
  
  if(hitcoords[1] %in% feature_start:feature_stop){
    if(direction == "forward"){
      print('1')
      overlap = feature_stop - hitcoords[1]
      if(overlap > overlap_for){
        overlap_for <- overlap
      }
    }
    
    if(direction == "backward"){
      print('2')
      overlap = hitcoords[1] - feature_start
      if(overlap > overlap_for){
        overlap_for <- overlap
      }
    }
    
  }
  
  if(hitcoords[2] %in% feature_start:feature_stop){
    print('3')
    if(direction == "forward"){
      overlap = hitcoords[2] - feature_start
      if(overlap > overlap_rev){
        overlap_rev <- overlap
      }
    }
    if(direction == "backward"){
      print('4')
      overlap = feature_stop - hitcoords[2]
      if(overlap > overlap_rev){
        overlap_rev <- overlap
      }
    }
  }
  
  return(c(overlap_for, overlap_rev))
}

feature_start = 43247
feature_stop = 44608
hitcoords = c(43629, 42987)
overlaps = c(0, 0)
feat_overlap(feature_start, feature_stop, hitcoords, direction = "backward", overlaps)


remove_dashes <- function(query){
  new_query <- query
  for(i in 1:length(new_query)){
    if(substr(new_query[i], 1,1) == '-'){
      new_query[i] <- str_replace_all(as.character(query[i]), '-', '')
    }
  }
  return(new_query)
}
remove_dashes(alignment[[2]])

#Function to blast the consensus sequence against N222.2 and obtain the coordinates of overlap, if any
#One issue is this isnt checking for cases where there is an entire gene contained within the igr. extremely unlikely but possible, and would be missed
#Also assumes overlaps are contiguous from either the start or the end. 
blast_check <- function(query, locus, blast_db){
  query <- remove_dashes(query)
  hits <- predict(blast_db, query)
  print(hits)
  hitcoords <- c(hits$sstart[1], hits$send[1])
  overlaps <- c(0, 0)
  direction <- "forward"
  
  if(hitcoords[1] > hitcoords[2]){
    direction <- "backward"
  }
  for(i in 2:nrow(reference_anno)){ #skip first row as that feature is the whole file
    feature_start <- reference_anno$start[i]
    feature_stop <-reference_anno$end[i]
    
    overlaps <- feat_overlap(feature_start, feature_stop, hitcoords, direction, overlaps)

  }
  
  return(overlaps)
}

blast_check(query = alignment[[2]], locus = "NEIS0320", blast_db = blast_db)
feat_overlap()
alignment[[2]][9] <- paste0('GGGGGGGG', alignment[[2]][9])
alignment

#Function to substr igrs and check if the shortened sequences still contain variation
#assumes all igrs given are the same length. need to make sure already shorter igrs arent shortened. 
#solution = substr() the alignment rather than the igrs themselves. then just see if the strings equal each other. 
igr_shorten <- function(alignment, overlaps, locus){
  igrs_new <- substr(alignment, (overlaps[1]+1), (nchar(alignment)-overlaps[2]))
  #check if still variable
  current_car <- hits[which(hits$Locus == locus), "Group"]
  
  newcar <- current_car
  #CAREFUL OF NA'S HERE
  uniques <- unique(igrs_new)
  if(length(uniques) == 1){
     newcar <- "NoVar"
  }
  if(length(uniques) == 2){
    if(length(alignment[which(alignment == uniques[1])]) >= 1 && length(alignment[which(alignment == uniques[2])]) == 1 && names(alignment[which(alignment == uniques[2])]) == paste0(locus, "_28269")){
      newcar <- "N188Var"
    }
    if(length(alignment[which(alignment == uniques[2])]) >= 1 && length(alignment[which(alignment == uniques[1])]) == 1 && names(alignment[which(alignment == uniques[1])]) == paste0(locus, "_28269")){
      newcar <- "N188Var"
    }
    else{
      newcar <- "Var"
    }
  }
  if(length(uniques) >2){
    newcar <- "Var"
  }
  
  change <- TRUE
  if(current_car == newcar){
    change <- FALSE
  }
  print(newcar)
  print(length(uniques))
  return(list(igrs_new, newcar, change))
}
overlaps <- c(0, 0)
igr_shorten(alignment[3:10], overlaps_up, locus = "NEIS0320")
overlaps_up


#Function to check if previous category of that locus has changed (i.e. used to be var or N188only but has changed to N188only or novar).
writer <- function(shortened_up, shortened_down, newhits, locus){
  cat_up <- shortened_up[2]
  cat_down <- shortened_down[2]
  print(cat_up)
  print(cat_down)
  change_up <- shortened_up[3]
  change_down <- shortened_down[3]
  change <- FALSE
  cat <- ''
  
  if(change_up == TRUE | change_down == TRUE){
    change <- TRUE
  }
  
  if(cat_up == 'NoVar' && cat_down == 'NoVar'){
    cat <- cat_up
  }
  if(cat_up == "N188Var" | cat_down == "N188Var"){
    cat <- "N188Var"
  }
  if(cat_up == "Var" | cat_down == "Var"){
    cat <- "Var"
  }
  print(paste("cat_up is:", cat_up))
  print(paste("cat_down is:", cat_down))
  print(cat)
  newhits[which(newhits$Locus == locus), "Change"] <- change
  newhits[which(newhits$Locus == locus), "Group" ] <- cat
  
  return(newhits)
}

main <- function(){
  #Set the N222.2.fasta file as a blast database
  up_igrs <- read.fasta("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/final_igrs_up.xmfa", as.string = TRUE)
  down_igrs <- read.fasta("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/final_igrs_down.xmfa", as.string = TRUE)
  hits <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/significant_loci_matches.xlsx")[,c(1:4,6,7)]
  blast_db <- blast(db="/Users/u5501917/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/N222_database/N222.2.fasta")
  reference_anno <- read_gbk("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/N222.2.gbk")
  reference_fas <- read.fasta("~/Documents/PhD/PhD/RNA_IGR/Isolate_fastas/N222.2.fasta", as.string = TRUE)
  newhits <- hits
  newhits$Change <- FALSE
  
  for(i in 1:nrow(hits)){
    i = 1
    locus <- hits$Locus[i]
    if(hits$Group[i] == "NoVar"){
      next
    }
    print(locus)
    #first extract igrs for a locus
    igrs <- extractor(locus = locus)
    #align igrs
    alignment <- align(igrs, locus)
    print("alignment successful")
    #identify region of overlap with coding sequences, if any
    overlaps_up <- blast_check(query = alignment[[1]], locus = locus, blast_db = blast_db)
    overlaps_down <- blast_check(query = alignment[[2]], locus = locus, blast_db = blast_db)
    print("overlaps identified")
    #Shorten igrs by the overlap and check if their category changes
    new_up <- igr_shorten(alignment = alignment[[1]][-1], overlaps = overlaps_up, locus)
    new_down <- igr_shorten(alignment = alignment[[2]][-1], overlaps = overlaps_down, locus)
    print("igrs shortened")
    #add columns to hits, save as new file
    newhits <- writer(new_up, new_down, newhits, locus)
    print('data added')
    print("==========")
    print("==========")
  }
  
  write.xlsx(newhits, "~/Documents/PhD/PhD/RNA_IGR/significant_loci_matches_postprocess2.xlsx")
}

main()


