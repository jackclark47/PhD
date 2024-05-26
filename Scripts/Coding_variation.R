#Script to load in Genome comparator data and check if genes with coding region variation also have FR variation.

library(openxlsx)

#load coding variable genes from GC data
coding_var <- read.xlsx("~/Downloads/all8_all_loci_genome_comparator.xlsx", sheet = 2, colNames = TRUE)
#2217 loci total
nrow(coding_var)/2217 * 100 #38.0% of all loci have coding region variation


#Load FR variable genes from FR data
FR_var <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx", sheet = 3)

#simple loop to check presence of coding_var loci in FR_var 
present_count = 0
present <- c()
absent <- c()
for(i in 1:nrow(coding_var)){
  if(coding_var$Locus[i] %in% FR_var$Locus){
    present_count = present_count + 1
    present <- c(present, coding_var$Locus[i])
  }
  else{
    absent <- c(absent, coding_var$Locus[i])
  }
}
present_count

present_count/nrow(FR_var) * 100
#56.7% of FR variable loci also have coding region variation
present_count/nrow(coding_var) * 100
#75.6% of coding region variable loci also have FR variation

#So coding region variation tends to be accompanied by FR variation but FR variation is substantially more common


#Now look into this in the context of RNAseq data - do significant DEGs tend to have coding variation?

#Load RNAseq data
rnaseq <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx", sheet = 8)
#What proportion of significant hits have coding region variation?
sigrna <- rnaseq[which(rnaseq$`log2(EXP/Min)` >= 1),]
sigrna <- sigrna[which(sigrna$X13 != 1),]

#what proportion of coding variable loci are significant?
count = 0
for(i in 1:nrow(sigrna)){
  if(sigrna$`1717.genes`[i] %in% coding_var$Locus){
    count = count + 1
  }
}
count

count/nrow(sigrna) *100 #29.2% of rnaseq loci that are significant have coding variation
count/nrow(coding_var) *100 #13.6% of coding variable loci are significant. 
nrow(sigrna)/2217 * 100 #17.6% of all loci are significant degs

#Do coding region variable loci have a significant different in LFC compared to non coding region variable loci? what about to FR variable loci?
  #do loci with both variations and compare against loci only with coding and loci only with FR variation



841+1300

1131/2217*100
