###Script to make pairwise distribution matrices for both RNAseq and IGR data
#Talk to Jonathan about this and ask if your outline makes sense

#load packages
library(openxlsx)

#Load data
setwd("~/Documents/PhD/PhD/RNA_IGR")
rna <- read.xlsx("11-7_run_5_transcripts-cdb-16082022.xlsx", sheet = 10)
igr <- read.xlsx("NEIS_IGR_Alleles-cdb-21092022.xlsx", sheet = 4)

###Matrix structure
#Each row will be a locus, each column a pairwise comparison between two isolates
#Must make sure both matrices have the same ordering
#Value will be 0 or 1 
#            N222.2_N457.1     N222.2_N455.3
#NEIS0001         1                   0
#NEIS0002         1                   1

#For IGR matrix, 1s and 0s represent same igr allele or different allele
#For RNA matrix, 1s and 0s represent same expression or different expression (maybe 0 is no significant log2fold change, 1 is significant)

###IGR matrix





###RNAseq matrix





###Pairwise comparisons between matrices
#two way ANOVA testing on each cell, colour significant cells green on a spreadsheet.


write.xlsx(x, "distributionmatrices.xlsx") #Mmake sheet of matrix coloured by significance as judged by ANOVA testing.