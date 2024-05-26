###Script to run linear regression on the igr vs rna datasets to find significant results

# the linear regression equation is:
# Y = a+bX+c + alpha
# 
# Y = vector of RNA transcription. 8 elements, one per isolate for a locus
# X = genotypes in all IGRs, coded as 0, 1, where 0 is the reference allele

#Load libraries
library(ggplot2)
library(openxlsx)
library(car)
library(stringr)
#Load in data
rna <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx', sheet = 6)[, c(8,9,15:22, 25:32)]
up_igr <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 1)[,c(1,3, 12:19,21)]
down_igr <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 1)[,c(1,2,4:11, 20)]

#Remove whitespaces in locus names in the rna dataset
for(i in 1:length(rna$`1717.genes`)){
  rna$`1717.genes`[i] <- str_trim(rna$`1717.genes`[i])
}

sig_up <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 3)[,c(1,3, 12:19,21)]
sig_down <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 3)[,c(1,2, 4:11,20)]

count = 0
for(i in 1:length(rna$`1717.genes`)){
  if(rna$`1717.genes`[i] %in% sig_up$Locus){
    count = count + 1
  }
}
count


#Order of isolates in both datasets:
#RNA: 20578 22698 N188.1  N222.1, N222.2, N445.1, N459.3, N459.6
#IGR: 20578 22698 N222.1  N188.1  N445.1, N222.2  N459.3, N459.6

#So rearrange columns in the rna dataset so theyre the same order
rna <- rna[,c(1:4, 6,5, 8,7, 9:12, 14,13,16,15, 17, 18)]

#For each locus build a table of:

# Isolate   expression    Allele_up   Allele_Down


#get an individual test case
isolates <- c(27509, 27553, 28262, 28269, 28287, 53930, 53948, 53951)

#Function that returns a table for a particular locus of the format described above. 
tabulator <- function(locus_name){
  if(!locus_name %in% up_igr$Locus){
    return(NA)
  }
    
  out <- as.data.frame(matrix(data= NA, nrow = 8, ncol = 4))
  colnames(out) <- c("Isolate", "log2_expression", "Allele_Up", "Allele_Down")
  out$Isolate <- isolates
  out$log2_expression <- as.vector(unlist(rna[which(rna$`1717.genes`==locus_name), 11:18][1,]))
  out$Allele_Up <- as.vector(unlist(up_igr[which(up_igr$Locus==locus_name), 3:10]))
  out$Allele_Down <- as.vector(unlist(down_igr[which(down_igr$Locus==locus_name), 3:10]))
  
  return(out)
}
tabulator("NEIS1769")
tabulator("NEIS1946")
#Create the empty dataframe to store our p values in
sigtable <- as.data.frame(matrix(data= NA, nrow = nrow(rna), ncol = 3))
colnames(sigtable) <- c('Locus', 'p_up', 'p_down')

count = 0
#A loop to populate the dataframe
for(i in 2:length(rna$`1717.genes`)){
  print(i)
  #get the name of the current locus and run the tabulator function
  locus_name <- str_trim(rna$`1717.genes`[i]) #remove emptyspaces in the rna dataset that sometimes appear
  out <- tabulator(locus_name)
  
  #enter the locus name into the empty df
  sigtable$Locus[i] <- locus_name
  
  #If the tabulator function returns an empty or shortened table, skip this locus.
  if(length(out) == 1){
    if(is.na(out)){
      count = count + 1
      next
    }
  }
  #if table contains at least one NA value in the expression columns, skip the locus
  if(sum(is.na(out$log2_expression)) > 0){
    count = count + 1
    next
  }

  #if the table ihas no NAs and has multiple alleles in its up igr, run a linear regression model on the locus
  if(length(unique(out$Allele_Up)) > 1){
    
    pval_up <- summary(lm(log2_expression ~ Allele_Up, data = out, na.action = na.omit))$coefficients
    if(length(pval_up != 8)){
      sigtable$p_up[i] <- NA
    }
    if(length(pval_up) == 8){
      print('yay')
      pval_up <- pval_up[2,4]
      sigtable$p_up[i] <- pval_up
    }
  }
  #likewise for the down_igr
  if(length(unique(out$Allele_Down)) > 1){
    
    pval_down <- summary(lm(log2_expression ~ Allele_Down, data = out, na.action = na.omit))$coefficients
    if(length(pval_down != 8)){
      sigtable$p_down[i] <- NA
    }
    
    if(length(pval_down) == 8){
      pval_down <- pval_down[2,4]
      print(pval_down)
      sigtable$p_down[i] <- pval_down
      print(sigtable$p_down[i])
    }

  }
  else{
    count = count + 1
  }
  
  #Reset values of p_down and p_up before the next loop
  pval_down <- 0
  pval_up <- 0
  
}
count
2272-1635 #637
#order the significance table by those loci which are significant in the up_igr comparison
sigtable <- sigtable[order(sigtable$p_up),]
isolates <- c(27509, 27553, 28262, 28269, 28287, 53930, 53948, 53951)

#If both p_up and p_down are NA values, remove them from the table
removables <- c()
for(i in 1:length(sigtable$Locus)){
  if(is.na(sigtable$p_up[i]) && is.na(sigtable$p_down[i])){
    removables <- c(removables, i)
  }
}

sigtable <- sigtable[-removables,]


#add igr and expression data to sigtable
for(i in 1:length(sigtable$Locus)){
  locusname <- sigtable$Locus[i]
  print(locusname)
  sigtable$upigr_27509[i] <- up_igr$`27509_up`[which(up_igr$Locus==locusname)]
  sigtable$upigr_27553[i] <- up_igr$`27553_up`[which(up_igr$Locus==locusname)]
  sigtable$upigr_28262[i] <- up_igr$`28262_up`[which(up_igr$Locus==locusname)]
  sigtable$upigr_28269[i] <- up_igr$`28269_up`[which(up_igr$Locus==locusname)]
  sigtable$upigr_28287[i] <- up_igr$`28287_up`[which(up_igr$Locus==locusname)]
  sigtable$upigr_53930[i] <- up_igr$`53930_up`[which(up_igr$Locus==locusname)]
  sigtable$upigr_53948[i] <- up_igr$`53948_up`[which(up_igr$Locus==locusname)]
  sigtable$upigr_53951[i] <- up_igr$`53951_up`[which(up_igr$Locus==locusname)]
  sigtable$downigr_27509[i] <- down_igr$`27509_down`[which(down_igr$Locus==locusname)]
  sigtable$downigr_27553[i] <- down_igr$`27553_down`[which(down_igr$Locus==locusname)]
  sigtable$downigr_28262[i] <- down_igr$`28262_down`[which(down_igr$Locus==locusname)]
  sigtable$downigr_28269[i] <- down_igr$`28269_down`[which(down_igr$Locus==locusname)]
  sigtable$downigr_28287[i] <- down_igr$`28287_down`[which(down_igr$Locus==locusname)]
  sigtable$downigr_53930[i] <- down_igr$`53930_down`[which(down_igr$Locus==locusname)]
  sigtable$downigr_53948[i] <- down_igr$`53948_down`[which(down_igr$Locus==locusname)]
  sigtable$downigr_53951[i] <- down_igr$`53951_down`[which(down_igr$Locus==locusname)]
  sigtable$rna_27509[i] <- rna$`log2(EXP/MIN)`[which(down_igr$Locus==locusname)]
  sigtable$rna_27553[i] <- rna$`X26`[which(down_igr$Locus==locusname)]
  sigtable$rna_28262[i] <- rna$`X28`[which(down_igr$Locus==locusname)]
  sigtable$rna_28269[i] <- rna$`X27`[which(down_igr$Locus==locusname)]
  sigtable$rna_28287[i] <- rna$`X30`[which(down_igr$Locus==locusname)]
  sigtable$rna_53930[i] <- rna$`X29`[which(down_igr$Locus==locusname)]
  sigtable$rna_53948[i] <- rna$`X31`[which(down_igr$Locus==locusname)]
  sigtable$rna_53951[i] <- rna$`X32`[which(down_igr$Locus==locusname)]
}

#Write to excel file
write.xlsx(sigtable, 'linearregression_igrRNA_pvalues_test.xlsx')

length(sigtable$p_up[which(sigtable$p_up <= 0.05)])
length(sigtable$p_down[which(sigtable$p_down <= 0.05)])

length(sigtable$Locus[which(sigtable$p_up <= 0.05 | sigtable$p_down <= 0.05)])

###Redo with 28269 excluded

rna <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx', sheet = 6)[, c(8,9,15:22, 25:32)]
up_igr <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 1)[,c(1,3, 12:19,21)]
down_igr <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 1)[,c(1,2,4:11, 20)]

#Remove whitespaces in locus names in the rna dataset
for(i in 1:length(rna$`1717.genes`)){
  rna$`1717.genes`[i] <- str_trim(rna$`1717.genes`[i])
}

#Order of isolates in both datasets:
#RNA: 20578 22698 N188.1  N222.1, N222.2, N445.1, N459.3, N459.6
#IGR: 20578 22698 N222.1  N188.1  N445.1, N222.2  N459.3, N459.6

#So rearrange columns in the rna dataset so theyre the same order
rna <- rna[,c(1:4, 6,5, 8,7, 9:12, 14,13,16,15, 17, 18)]

#Remove 28269 aka N188.1
rna <- rna[,-c(6)]

#For each locus build a table of:

# Isolate   expression    Allele_up   Allele_Down


#get an individual test case
isolates <- c(27509, 27553, 28262, 28287, 53930, 53948, 53951)

#Function that returns a table for a particular locus of the format described above. 
tabulator <- function(locus_name){
  if(!locus_name %in% up_igr$Locus){
    return(NA)
  }
  
  out <- as.data.frame(matrix(data= NA, nrow = 7, ncol = 4))
  colnames(out) <- c("Isolate", "log2_expression", "Allele_Up", "Allele_Down")
  out$Isolate <- isolates
  out$log2_expression <- as.vector(unlist(rna[which(rna$`1717.genes`==locus_name), 11:17][1,]))
  out$Allele_Up <- as.vector(unlist(up_igr[which(up_igr$Locus==locus_name), 3:9]))
  out$Allele_Down <- as.vector(unlist(down_igr[which(down_igr$Locus==locus_name), 3:9]))
  
  return(out)
}
tabulator("NEIS1769")
tabulator("NEIS1967")
#Create the empty dataframe to store our p values in
sigtable <- as.data.frame(matrix(data= NA, nrow = nrow(rna), ncol = 3))
colnames(sigtable) <- c('Locus', 'p_up', 'p_down')

#A loop to populate the dataframe
for(i in 2:length(rna$`1717.genes`)){
  print(i)
  #get the name of the current locus and run the tabulator function
  locus_name <- str_trim(rna$`1717.genes`[i]) #remove emptyspaces in the rna dataset that sometimes appear
  out <- tabulator(locus_name)
  
  #enter the locus name into the empty df
  sigtable$Locus[i] <- locus_name
  
  #If the tabulator function returns an empty or shortened table, skip this locus.
  if(length(out) == 1){
    if(is.na(out)){
      next
    }
  }
  #if table contains at least one NA value in the expression columns, skip the locus
  if(sum(is.na(out$log2_expression)) > 0){
    next
  }
  
  
  if(length(unique(out$Allele_Up)) > 1){
    
    pval_up <- summary(lm(log2_expression ~ Allele_Up, data = out, na.action = na.omit))$coefficients
    if(length(pval_up != 8)){
      sigtable$p_up[i] <- NA
    }
    if(length(pval_up) == 8){
      print('yay')
      pval_up <- pval_up[2,4]
      sigtable$p_up[i] <- pval_up
    }
  }
  
  
  
  #if the table ihas no NAs and has multiple alleles in its up igr, run a linear regression model on the locus
  if(length(unique(out$Allele_Up)) > 1){
    pval_up <- summary(lm(log2_expression ~ Allele_Up, data = out, na.action = na.omit))$coefficients
    
    if(length(pval_up != 8)){
      sigtable$p_up[i] <- NA
    }
    if(length(pval_up) == 8){
      print('yay')
      pval_up <- pval_up[2,4]
      sigtable$p_up[i] <- pval_up
    }
  }
  
  #likewise for the down_igr
  if(length(unique(out$Allele_Down)) > 1){
    pval_down <- summary(lm(log2_expression ~ Allele_Down, data = out, na.action = na.omit))$coefficients
    
    if(length(pval_down != 8)){
      sigtable$p_down[i] <- NA
    }
    if(length(pval_down) == 8){
      print('yay')
      pval_down <- pval_down[2,4]
      sigtable$p_down[i] <- pval_down
    }
  }
  
  #Reset values of p_down and p_up before the next loop
  pval_down <- 0
  pval_up <- 0
  
}

#order the significance table by those loci which are significant in the up_igr comparison
sigtable <- sigtable[order(sigtable$p_up),]
isolates <- c(27509, 27553, 28262, 28287, 53930, 53948, 53951)

#If both p_up and p_down are NA values, remove them from the table
removables <- c()
for(i in 1:length(sigtable$Locus)){
  if(is.na(sigtable$p_up[i]) && is.na(sigtable$p_down[i])){
    removables <- c(removables, i)
  }
}

sigtable <- sigtable[-removables,]

#add igr and expression data to sigtable
for(i in 1:length(sigtable$Locus)){
  locusname <- sigtable$Locus[i]
  print(locusname)
  sigtable$upigr_27509[i] <- up_igr$`27509_up`[which(up_igr$Locus==locusname)]
  sigtable$upigr_27553[i] <- up_igr$`27553_up`[which(up_igr$Locus==locusname)]
  sigtable$upigr_28262[i] <- up_igr$`28262_up`[which(up_igr$Locus==locusname)]
  sigtable$upigr_28287[i] <- up_igr$`28287_up`[which(up_igr$Locus==locusname)]
  sigtable$upigr_53930[i] <- up_igr$`53930_up`[which(up_igr$Locus==locusname)]
  sigtable$upigr_53948[i] <- up_igr$`53948_up`[which(up_igr$Locus==locusname)]
  sigtable$upigr_53951[i] <- up_igr$`53951_up`[which(up_igr$Locus==locusname)]
  sigtable$downigr_27509[i] <- down_igr$`27509_down`[which(down_igr$Locus==locusname)]
  sigtable$downigr_27553[i] <- down_igr$`27553_down`[which(down_igr$Locus==locusname)]
  sigtable$downigr_28262[i] <- down_igr$`28262_down`[which(down_igr$Locus==locusname)]
  sigtable$downigr_28287[i] <- down_igr$`28287_down`[which(down_igr$Locus==locusname)]
  sigtable$downigr_53930[i] <- down_igr$`53930_down`[which(down_igr$Locus==locusname)]
  sigtable$downigr_53948[i] <- down_igr$`53948_down`[which(down_igr$Locus==locusname)]
  sigtable$downigr_53951[i] <- down_igr$`53951_down`[which(down_igr$Locus==locusname)]
  sigtable$rna_27509[i] <- rna$`log2(EXP/MIN)`[which(down_igr$Locus==locusname)]
  sigtable$rna_27553[i] <- rna$`X26`[which(down_igr$Locus==locusname)]
  sigtable$rna_28262[i] <- rna$`X28`[which(down_igr$Locus==locusname)]
  sigtable$rna_28287[i] <- rna$`X30`[which(down_igr$Locus==locusname)]
  sigtable$rna_53930[i] <- rna$`X29`[which(down_igr$Locus==locusname)]
  sigtable$rna_53948[i] <- rna$`X31`[which(down_igr$Locus==locusname)]
  sigtable$rna_53951[i] <- rna$`X32`[which(down_igr$Locus==locusname)]
}

sigtable2 <- sigtable
sigtable2$test <- c(1)
#Write to excel file
write.xlsx(sigtable, 'linearregression_igrRNA_pvalues_no28269_test.xlsx')

length(sigtable$p_up[which(sigtable$p_up <= 0.05)])
length(sigtable$p_down[which(sigtable$p_down <= 0.05)])

length(sigtable$Locus[which(sigtable$p_up <= 0.05 | sigtable$p_down <= 0.05)])

