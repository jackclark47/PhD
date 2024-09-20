#Script to assess linear regression outputs. Sees how many linreg hits are also hits in the rna dataset
library(openxlsx)

linreg <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/linearregression_igrRNA_pvalues_test.xlsx")
#linreg <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/linearregression_igrRNA_pvalues_no28269_test.xlsx")
rna <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx", sheet = 8)
#Remove whitespaces in locus names in the rna dataset
for(i in 1:length(rna$`1717.genes`)){
  rna$`1717.genes`[i] <- str_trim(rna$`1717.genes`[i])
}


#How many linreg hits are in the rna dataset?
count = 0
loci <- c()
for(i in linreg$Locus){
  p_up = linreg[which(linreg$Locus == i),2][1]
  p_down = linreg[which(linreg$Locus == i),3][1]
  if(is.na(p_up)){
    p_up = 100
  }
  if(is.na(p_down)){
    p_down = 100
  }
  if(i %in% rna$`1717.genes` && (p_up < 0.05 | p_down < 0.05)){
    count = count + 1
    loci <- c(loci, i)
  }
}
#142 significant linreg hits are in the rna dataset overall
#How many of these have logfoldchange >1 for MAX/MIN?
linreg_hits <- linreg[which(linreg$Locus %in% loci),]
count = 0
loci <- c()
for(i in linreg_hits$Locus){
  logchange <- rna[which(rna$`1717.genes` == i),14]
  if(logchange >= 1){
    count = count + 1
    loci <- c(loci, i)
  }
}
#103 (72.5%) linreg hits have logfold change MAX/MIN >= 1

#And how many are significant in the rna?
linreg_changed <- linreg_hits[which(linreg_hits$Locus %in% loci),]
count = 0
loci <- c()
for(i in linreg_changed$Locus){
  qval <- rna[which(rna$`1717.genes` == i),13]
  if(qval < 0.05){
    count = count + 1
    loci <- c(loci, i)
  }
}
#34 (23.9%)linreg hits are significantly differently expressed in the rna dataset when using MAX/MIN



#How many linear regression hits are signficantly differently expressed in any pairwise comparison?

rna[which(rna$`1717.genes` == 'NEIS1594'),]
for(i in rna[1359,101:156]){
  if(i < 0.05){
    count = count + 1
    loci <- c(loci, i)
    break
  }
}

count = 0
loci <- c()
for(i in linreg_changed$Locus){
  for(j in rna[which(rna$`1717.genes` == i),101:156]){
    if(j < 0.05){
      count = count + 1
      loci <- c(loci, i)
      break
    }
  }
}

hits <- linreg[which(linreg$Locus %in% loci),]
hits <- hits[,c(1,5:8)]
