#script to concatenate all igrs for an isolate, performed for each isolate
###Up and Down IGRs are stored separately, so we need to alternate between them

#Libraries
library(seqinr)
library(msa)
###FIRST NEED TO MAKE SURE IM ONLY LOOKING AT IGRS PRESENT IN EVERY ISOLATE

#Remove loci which contain at least one NA / have a size of 0 


#Make function
seq_join <- function(start){
  
  #Up_IGRs
  up_igrs <- read.fasta("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/final_igrs_up.xmfa", as.string = TRUE)[seq(start, 15568, 8)]
  #Down_IGRs
  down_igrs <- read.fasta("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/final_igrs_down.xmfa", as.string = TRUE)[seq(start, 15568, 8)]

  
  full_seq <- ""
  for(i in 1:(length(up_igrs))){
    full_seq <- paste(full_seq, up_igrs[i], sep = "")
    full_seq <- paste(full_seq, down_igrs[i], sep = "")
    
  }
  return(full_seq)
}

#Run function
i_27509 <- seq_join(1)
i_27553 <- seq_join(2)
i_28262 <- seq_join(3)
i_28269 <- seq_join(4)
i_28287 <- seq_join(5)
i_53930 <- seq_join(6)
i_53948 <- seq_join(7)
i_53951 <- seq_join(8)

#combined in a single list
seqs <- list(i_27509, i_27553, i_28262, i_28269, i_28287, i_53930, i_53948, i_53951)
seqnames <- c("27509", "27553", "28262", "28269", "28287", "53930", "53948", "53951")
names(seqs) <- seqnames
#this writes a file containing 8 sequences, each of which is the convatenated core up and downstream igrs for a given isolate. 
#File is written in locus order, with the up igr, then down igr, then next locus up_igr, down_igr etc
write.fasta(seqs, names = c("27509", "27553", "28262", "28269", "28287", "53930", "53948", "53951"), "~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igrs/Concatenated_core_igrs/all_seqs_igrs.fasta")
 

#Try aligning with Clustal within R:
msa(seqs, method="ClustalOmega", type="dna")



