#Main function, first converts a DNAStringSet to a vector of sequences using DNAStringtoCharacter(),
#then creates a list of every pairwise combination of strings (a bit like every possible tuple),
#then calculates the number of character differences between every single pairwise combination using diversity()
#Adds those all together and divides the total by the number of pairwise combinations to get pi i.e. nucleotide diversity

my.nuc.div <- function(sequences){ #takes output from AlignSeqs as input, which is a DNAStringSet object
  sequences <- DNAStringtoCharacter(sequences) 
  p.diversities <- c()
  seq_pairs <- combn(sequences, 2)
  for(i in 1:length(seq_pairs[1,])){
    s1 <- unlist(seq_pairs[1,i])
    s2 <- unlist(seq_pairs[2,i])
    p.diversities <- c(p.diversities, diversity(s1,s2))
  }
  pi <- sum(p.diversities)/length(seq_pairs[1,])
  return(pi)
}

#Function that finds the number of differences between each element of a DNA sequence pair
diversity <- function(s1, s2){
  count = 0
  for(char in 1:nchar(s1)){
    if(substr(s1,char,char) != substr(s2,char,char)){
      count = count + 1
    }
  }
  diversity = count/nchar(s1) 
  return(diversity)
}

#Function that converts DNAString objects from ape to a vector of strings for easier handling
DNAStringtoCharacter <- function(DNAStringobject){ #Input is a DNAStringobject output from AlignSeqs
  sequence_list <- c()
  len <- width(DNAStringobject)[1]
  sequences <- unlist(DNAStringobject) 
  sequences <- as.character(sequences)
  
  iterations <- nchar(sequences)/len
  for(i in 1:iterations){ 
    stop = len*i
    start = len*i - len+1
    sequence <- substr(sequences, start,stop)
    sequence_list <- c(sequence_list, sequence)
  }
  
  return(sequence_list)
}
