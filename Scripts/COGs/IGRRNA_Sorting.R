#Script to obtain lists of loci and conduct GO Term analysis

#Questions: 
#Am I just looking at significant loci in these groups? 
#Or all loci, so I can see a ratio of sig:nonsig for each group?

#Groups:
#HighVar, HighExp
#HighVar, LowExp
#LowVar, HighExp
#LowVar, LowExp

#load datasets

rna <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx', sheet = 6)[, c(8,9,15:22, 25:32)]
up_igr <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_12102023.xlsx', sheet = 1)[,c(1,3, 12:19,21)]
down_igr <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_12102023.xlsx', sheet = 1)[,c(1,2,4:11, 20)]
pi_values <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/pi_values.xlsx")
#Remove extra spaces in the rna locus names  
for(i in 1:length(rna$`1717.genes`)){
  rna$`1717.genes`[i] <- str_trim(rna$`1717.genes`[i])
}

#Definitions
#HighVar - Loci with more than x nucelotide diversity
#LowVar - Loci with less than or equal to x nucleotide diversity
#HighExp - Loci with a log2fold change greater than or equal to 1 in the EXP/MIN column in the rna dataset
#LowExp - Loci with an log2fold change lower than 1 in the EXP/MIN column in the rna dataset
  
###Identify loci belonging to each group
#Initialise variables
hi_exp <- c()
low_exp <- c()
for(i in 2:length(rna$`1717.genes`)){
  print(i)
  if(is.na(rna$`1717.genes`[i]) | is.na(rna$`log2(EXP/MIN)`[i])){
    next
  }
  if(rna$`log2(EXP/MIN)`[i] >= 1){
    hi_exp <- c(hi_exp, rna$`1717.genes`[i])
  }
  if(rna$`log2(EXP/MIN)`[i] < 1){
    low_exp <- c(low_exp, rna$`1717.genes`[i])
  }
}


#Change the below to take into account nucleotide diversity
hi_var <- c()
low_var <- c()
for(i in 1:length(up_igr$`Alleles-Up`)){
  print(i)
  if(is.na(up_igr$Locus[i]) | is.na(down_igr$`Alleles-Down`[i]) | is.na(up_igr$`Alleles-Up`[i])){
    next
  }
  if(sum(up_igr$`Alleles-Up`[i], down_igr$`Alleles-Down`[i]) > 4){
    hi_var <- c(hi_var, up_igr$Locus[i])
  }
  if(sum(up_igr$`Alleles-Up`[i], down_igr$`Alleles-Down`[i]) <= 4){
    low_var <- c(low_var, up_igr$Locus[i])
  }
}

#look for overlap across groups to define the final groupings. 
hivar_hiexp <- hi_var[which(hi_var %in% hi_exp)] #27 loci
hivar_lowexp <- hi_var[which(hi_var %in% low_exp)] #51 loci
lowvar_hiexp <- low_var[which(low_var %in% hi_exp)] #341 loci
lowvar_lowexp <- low_var[which(low_var %in% low_exp)] #1075 loci

#What proportion of loci in each group are significant in the linear regression model?
sigtable <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/linearregression_igrRNA_pvalues.xlsx")
hivar_hiexp[which(hivar_hiexp %in% sigtable$Locus)] #23.  so 85%
hivar_lowexp[which(hivar_lowexp %in% sigtable$Locus)] #47. so 92%
lowvar_hiexp[which(lowvar_hiexp %in% sigtable$Locus)] #80. so 24%
lowvar_lowexp[which(lowvar_lowexp %in% sigtable$Locus)]#149  so 14%


