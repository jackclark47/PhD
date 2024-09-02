#Script to get a list of known regulators in N meningitidis from existing data
#then idnetify which are contained in the rnaseq data
#Then read about downstream effects of the regulators and extract any of their downstream regulated loci

library(openxlsx)

regulator_init <- read.xlsx('~/Documents/PhD/PhD/Sequencing_and_Annotating/gene_classes.xlsx', sheet = 'Regulation')[-c(1,2),c(3:6)]
rna <- read.xlsx('~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx', sheet = 8)
fr <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 1)


#Fix rna NEIS code names
for(i in 1:length(rna$`1717.genes`)){
  rna$`1717.genes`[i] <- str_trim(rna$`1717.genes`[i])
}


regulator_init <- regulator_init[rowSums(is.na(regulator_init)) != ncol(regulator_init),]
nrow(regulator_init)
colnames(regulator_init) <- c('pubmlst', 'NMB', 'KO', 'description')
#Check presence of each regulator in the RNAseq data
#if an NMB code is available, use that, otherwise use NEIS. If no match found for NMB, try NEIS too.
rna_regulators <- c()
rna_indices <- c()
for(i in 1:nrow(regulator_init)){
  if(!is.na(regulator_init$NMB[i])){
    if(regulator_init$NMB[i] %in% rna$`1359.genes`){
      rna_regulators <- c(rna_regulators, i)
      rna_index <- which(rna$`1359.genes` == regulator_init$NMB[i])
    } else if(regulator_init$pubmlst[i] %in% rna$`1717.genes`){
      if(is.na(regulator_init$pubmlst[i])){
        next
      }
      rna_regulators <- c(rna_regulators, i)
      rna_index <- which(rna$`1717.genes` == regulator_init$pubmlst[i])
    }
    rna_indices <- c(rna_indices, rna_index)
    next
  }
  
  if(!is.na(regulator_init$pubmlst[i])){
    if(regulator_init$pubmlst[i] %in% rna$`1717.genes`){
      rna_regulators <- c(rna_regulators, i)
      rna_index <- which(rna$`1717.genes` == regulator_init$pubmlst[i])
      rna_indices <- c(rna_indices, rna_index)
    }
  }
}

regulator_rna <- regulator_init[rna_regulators,]
regulator_rna_absent <- regulator_init[-rna_regulators,]

#Check presence in igr data - can only use NEIS codes here
fr_regulators <- c()
for(i in 1:nrow(regulator_init)){
  if(!is.na(regulator_init$pubmlst[i])){
    if(regulator_init$pubmlst[i] %in% fr$Locus){
      fr_regulators <- c(fr_regulators, i)
    }
  }
}

regulator_fr <- regulator_init[fr_regulators,]
regulator_fr_absent <- regulator_init[-fr_regulators,]


#how many are present in both igr and rnaseq datasets?
length(fr_regulators) #230
length(rna_regulators) #243
sum(fr_regulators %in% rna_regulators) #201

#which ones are absent in one of them?
#regulator_rna_nofr <- regulator_rna[which(!(regulator_rna$pubmlst %in% regulator_fr)),]

#get fr data for regulators
fr_data <- fr[which(fr$Locus %in% regulator_fr$pubmlst),]

#get rna data for regulators
is.na(regulator_rna$pubmlst) <- 0
regulator_rna[which(is.na(regulator_rna$pubmlst)),] <- 0
regulator_rna[which(is.na(regulator_rna$NMB)),] <- 0

regulator_rna$pubmlst
rna_data <- rna[rna_indices,]



table(unique(fr_data$Locus %in% rna_data$`1717.genes`)) #198 overlap here - 3 less than before, maybe due to non NEIS codes in frdata not being in the 1717 gene column

#In just the rna_data, find significantly different regulators
rna_sigdata <- rna_data[which(rna_data$`log2(EXP/Min)` >= 1 & as.numeric(rna_data$X13) < 0.05),] #64 sig regulators

#in fr_data, find fr variable regulators
fr_var <- fr_data[which(fr_data$Combined > 2),] #103 regulators are fr variable
table(rna_sigdata$`1717.genes` %in% fr_var$Locus) #21 are fr variable
table(rna_sigdata$`1717.genes` %in% fr_data$Locus) #47 are in the fr dataset
21/47*100 #45% are fr variable


class(rna_data$X13)
class(1)

#get list of regulator NEIS codes to get genic variation information from pubmlst
NEIS_ids <- rna_data$`1717.genes`
NEIS_ids <- unique(NEIS_ids)
write.csv(NEIS_ids, '~/Documents/PhD/PhD/Regulators/pubmlst_ids_regulators.csv', quote = F, row.names = F, col.names = F)


gc_var <- read.xlsx('~/Documents/PhD/PhD/Regulators/Genome_comparator_on_regulators_in_rnaseq.xlsx', sheet = 2)
gc_novar <- read.xlsx('~/Documents/PhD/PhD/Regulators/Genome_comparator_on_regulators_in_rnaseq.xlsx', sheet = 3)
#how many gc_var overlap with fr?
table(gc_var$Locus %in% fr_var$Locus) #49/73 are also fr variable

#how many regulators are gc_var or fr_var?
nrow(fr_var) + 73-49 #127 of the regulators in are gc or fr variable

#make full table of regulators

#Locus    RNAseq?   sig?    frvar?    gcvar?    
#NEIS001    1         1       Downvar      0
#NEIS002    0         0       Novar        1

#1 or 0 is present/absent or yes/no. frvar also has the vargroup the locus belongs to
df <- as.data.frame(matrix(NA, nrow = nrow(regulator_init), ncol = 5))
colnames(df) <- c('Locus', 'RNAseq', 'significant', 'frvar', 'gcvar')
df$Locus <- regulator_init$pubmlst

#if regulatorinitpubmlst is na, use nmb insteadf
for(i in 1:nrow(df)){
  if(is.na(df$Locus[i])){
    df$Locus[i] <- regulator_init$NMB[i]
  }
}

df$RNAseq <- 0
for(i in 1:nrow(df)){
  if(df$Locus[i] %in% rna$`1717.genes` | df$Locus[i] %in% rna$`1359.genes`){
    df$RNAseq[i] <- 1
  }
}
table(df$RNAseq)  

df$significant <- 0
for(i in 1:nrow(df)){
  if(df$RNAseq[i] == 1){
    rna_index <- which(rna$`1717.genes` == df$Locus[i])
    if(length(rna_index) == 0){
      rna_index <- which(rna$`1359.genes` == df$Locus[i])
    }
  }
  
  if(length(rna_index) > 1){
    print(df$Locus[i])
    print(rna_index)
    rna_index <- rna_index[1]
  }

  rna_row <- rna[rna_index,]
  if(as.numeric(rna_row$X13) < 0.05 & rna_row$`log2(EXP/Min)` >= 1){
    print('yes')
    df$significant[i] <- 1
  }
  
  if(df$RNAseq[i] == 0){
    df$significant[i] <- NA
  }
}

table(df$significant)  

NoVar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 2)
DownVar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 5)
UpVar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 6)
LowVar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 4)
HighVar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 7)


df$frvar <- NA
for(i in 1:nrow(df)){
  if(df$Locus[i] %in% NoVar$Locus){
    df$frvar[i] <- 'NoVar'
  }
  if(df$Locus[i] %in% DownVar$Locus){
    df$frvar[i] <- 'DownVar'
  }
  if(df$Locus[i] %in% UpVar$Locus){
    df$frvar[i] <- 'UpVar'
  }
  if(df$Locus[i] %in% LowVar$Locus){
    df$frvar[i] <- 'LowVar'
  }
  if(df$Locus[i] %in% HighVar$Locus){
    df$frvar[i] <- 'HighVar'
  }
}
  
table(df$frvar)
29+33+3+38 == 103

df$gcvar <- NA
for(i in 1:nrow(df)){
  if(df$Locus[i] %in% gc_var$Locus){
    df$gcvar[i] <- 1
  }
  if(df$Locus[i] %in% gc_novar$Locus){
    df$gcvar[i] <- 0
  }
}

table(df$gcvar)
#number of non gcvar and non fr var regulators
nrow(df[df$gcvar != 1 & (df$frvar == 'NoVar' | is.na(df$frvar)),]) #203
nrow(df) #325
203/325 #62%

#save spreadsheet
write.xlsx(df, '~/Documents/PhD/PhD/Regulators/regulator_variation_summary.xlsx')







  