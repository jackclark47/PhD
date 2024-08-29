#script to check if the isolate with FR variation is the one with lowest expression.
library(openxlsx)
#load data
fr <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 1)
rna <- read.xlsx('~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx', sheet = 8)

#remove whitespace in rna
rna$`1717.genes`<-str_trim(rna$`1717.genes`)
rna$X13 <- as.numeric(rna$X13)
rna$`log2(EXP/Min)` <- as.numeric(rna$`log2(EXP/Min)`)


fr <- fr[which(fr$Locus %in% rna$`1717.genes`),]


single_frs <- function(fr, colnum){
  single_variable <- c()
  for(i in 1:nrow(fr)){
    down_allele <- fr[i,colnum]
    up_allele <- fr[i,colnum+8]
    
    try(if((sum(na.omit(fr[i,4:11]) == down_allele) == 1 && down_allele == 2 && fr$`Alleles-Down`[i] == 2) &&
       sum(na.omit(fr[i, 12:19]) == up_allele) == 1 && up_allele == 2 && fr$`Alleles-Up`[i] == 2){
      
      single_variable <- c(i, single_variable)
    })
    
    try(if(sum(na.omit(fr[i,4:11]) == down_allele) == 1 && fr$`Alleles-Up`[i] == 1 && fr$`Alleles-Down`[i] == 2){
      single_variable <- c(i, single_variable)
    })
    
    try(if(sum(na.omit(fr[i,12:19]) == up_allele) == 1 && fr$`Alleles-Down`[i] == 1 && fr$`Alleles-Up`[i] == 2){
      single_variable <- c(i, single_variable)
    })
  }
  single_fr <- fr[single_variable,]
  return(single_fr)
}

test <- single_frs(fr, 7)

#531 fr variable only in 28269
length(test)

single_rnas <- function(single_fr, rna){
  single_rna <- rna[which(rna$`1717.genes` %in% single_fr$Locus),]
  single_rna_sig <- single_rna[which(single_rna$`log2(EXP/Min)` >= 1),]
  single_rna_sig <- single_rna_sig[which(single_rna_sig$X13 < 0.05),]
  
  return(list('rnasig' = single_rna_sig, 'rna' = single_rna))
}

#which of these loci were significant in rna?
rna_28269 <- rna[which(rna$`1717.genes` %in% fr_28269$Locus),] #all 534 are in rnaseq
rna_28269_sig <- rna_28269[which(rna_28269$`log2(EXP/Min)` >= 1),]
rna_28269_sig <- rna_28269_sig[which(rna_28269_sig$X13 < 0.05),] #112 are significant in rnaseq so about 65 are var in other loci

minmax_rnas <- function(single_rna, colnumrna){
  minmax <- c()
  for(i in 1:nrow(single_rna)){
    counts <- single_rna[i, colnumrna]
    min_exp <- single_rna$X23[i]
    max_exp <- max(as.numeric(single_rna[i, 15:22]))
    
    if(counts == min_exp | counts == max_exp){
      minmax <- c(minmax, i)
    }
  }
  rna_minmax <- single_rna[minmax,]
  return(rna_minmax)
}

#which of these loci were lowest/highest expressed in 28269?
minmax <- c()
for(i in 1:nrow(rna_28269)){
  counts_28269 <- rna_28269$X17[i]
  min_exp <- rna_28269$X23[i]
  max_exp <- max(as.numeric(rna_28269[i, 15:22]))

  if(counts_28269 == min_exp | counts_28269 == max_exp){
    minmax <- c(minmax, i)
  }
}

rna_28269_minmax <- rna_28269[minmax,]
nrow(rna_28269_minmax)
#of the 582 loci, 362 are min or max in 28269.

#match the two up
rna_crossover <- rna_28269_sig[which(rna_28269_sig$`1717.genes` %in% rna_28269_minmax$`1717.genes`),]
nrow(rna_crossover) #66 are sig and min/max expressed in 28269

list('test' = 1, 'other' = 3)

main <- function(colnum, colnumrna, fr, rna){
  
  single_fr <- single_frs(fr, colnum)
  single_variable_frs <- nrow(single_fr)
  print(single_variable_frs)
  
  single_rnas_out <- single_rnas(single_fr, rna)
  
  single_rna <- single_rnas_out$rna
  single_rna_sig <- single_rnas_out$rnasig
  
  num_in_rna <- nrow(single_rna)
  num_sig_in_rna <- nrow(single_rna_sig)
  print(num_in_rna)
  print(num_sig_in_rna)
  
  minmax_rna <- minmax_rnas(single_rna, colnumrna)
  num_minmax <- nrow(minmax_rna)
  print(num_minmax)
  
  crossover_rna <- single_rna_sig[which(single_rna_sig$`1717.genes` %in% minmax_rna$`1717.genes`),]
  num_crossover <- nrow(crossover_rna)
  print(num_crossover)
  
  category <- c('single_variable_frs', 'num_in_rna', 'num_sig_in_rna', 'num_minmax', 'num_crossover')
  value <- c(single_variable_frs, num_in_rna, num_sig_in_rna, num_minmax, num_crossover)
  
  out <- as.data.frame(cbind(category, value))

  return(out)
}


i27509 <- main(4,15, fr, rna) 
i27553 <- main(5,16, fr, rna)
i28262 <- main(6,18, fr, rna)
i28269 <- main(7,17, fr, rna)
i28287 <- main(8,20, fr, rna)
i53930 <- main(9,19, fr, rna)
i53948 <- main(10,21, fr, rna)
i53951 <- main(11,22, fr, rna)

out <- cbind(i27509, i27553$value, i28262$value, i28269$value, i28287$value, i53951$value)
out$`53930` <- 0
out$`53948` <- 0
colnames(out) <- c('category', '27509', '27553', '28262', '28269', '28287', '53951', '53930', '53948')
