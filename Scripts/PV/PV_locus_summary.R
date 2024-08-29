#function that makes a short summary table of PV and expression data for a given locus
library(openxlsx)
library(stringr)
#Load PV data
PV <- read.xlsx('~/Documents/PhD/PhD/PhasomeIt_data/phasomeit_out_RNAcheck.xlsx', sheet = 5)
#Load RNA data
Exp <- read.xlsx('~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx', sheet=8)[-1,]

#Preprocess
PV <- PV[which(!is.na(PV$PubMLST_id)),]
Exp$`1717.genes`<-str_trim(Exp$`1717.genes`)
Exp$X13 <- as.numeric(Exp$X13)
Exp$`log2(EXP/Min)` <- as.numeric(Exp$`log2(EXP/Min)`)

summarise <- function(locus, PV, rna){
  error = FALSE
  if(!(locus %in% rna$`1717.genes`)){
    warning(paste('Locus id ', locus, ' not in rnaseq data', sep = ''), immediate. = T)
    error = TRUE
  }
  if(!(locus %in% PV$PubMLST_id)){
    warning(paste('Locus id ', locus, ' not in PV data', sep = ''), immediate. = T)
    error = TRUE
  }
  
  if(error == TRUE){
    stop('Locus id ', locus, ' absent in rnaseq or PV data, or both', sep ='')
  }
  
  PVdata <- as.vector(unlist(as.vector(PV[which(PV$PubMLST_id == locus),20:27])))
  rnadata <- rna[which(rna$`1717.genes` == locus),15:22]
  rnadata <- as.vector(unlist(as.numeric(as.vector(rnadata[,c(1,2,4,3,6,5,7,8)]))))
  summary <- cbind(PVdata, rnadata)
  rownames(summary) <- c('27509', '27553', '28262', '28269', '28287', '53930', '53948', '53951')

  summary <- as.data.frame(summary)
  summary <- summary[order(summary$PVdata),]
  
  colnames(summary) <- c('pv_tract', 'expression')
  return(summary)
}

summarise('NEIS0568', PV, Exp)

