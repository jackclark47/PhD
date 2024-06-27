#PV analysis script
#Correlates PV state or IGR tract length with gene expression changes

#Libraries


#Load PV data
PV <- read.xlsx("~/Documents/PhD/PhD/PhasomeIt_data/phasomeit_out_RNAcheck.xlsx",sheet = 5)[1:54,c(1:4, 18:28,31:39)]

#Load expression data
Exp <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx",sheet = 8)[-c(1),]
#Remove whitespaces in locus names in the rna dataset
for(i in 1:length(Exp$`1717.genes`)){
  Exp$`1717.genes`[i] <- str_trim(Exp$`1717.genes`[i])
}

#Add NMB codes for opa genes to Exp
Exp[586,9] <- 'NMB0926'
Exp[1131,9] <- 'NMB1465'
Exp[1305,9] <- 'NMB1636'
Exp[1484,9] <- 'NMB0442'
Exp[1949,9] <- 'N222_01750'
Exp[230,9] <- 'N222_00206'

#Filter PV data to contain only rows where there is a difference in IGR tract length in carriage isolates, or a change in ON/OFF state for genic PV tracts
removables <- c()
IGRs <- c()
coding <- c()
for(i in 1:nrow(PV)){
  print(i)
  states <- unlist((PV[i,c(16:23)]))
  print(states)
  if(length(na.omit(states)) == 0){
    removables <- c(removables, i)
    next
  }
  if(unique(na.omit(states))[1] == 'IGR'){
    IGRs <- c(IGRs, i)
  }
  if(unique(na.omit(states))[1] != 'IGR'){
    coding <- c(coding, i)
  }
}
length(coding) #32 coding PV
length(IGRs) #21 IGR PV
length(removables) #1
nrow(PV) #54
#1 rows of 54 will be removed,leaving 53 PV loci (though some tracts are in the same locus)

#Split into IGR and coding PV tracts
PVfilt <- PV[-removables,]
PV_IGR <- PV[IGRs,]
PV_coding <- PV[coding,]

#for coding tracts, group by ON/OFF state, then test for significance with a wilcox test
out_coding <- as.data.frame(matrix(data = NA, nrow = nrow(PV_coding), ncol = 8))
colnames(out_coding) <- c('ID', 'Name', 'Function', 'PV_location', 'Status', 'Pval', 'MAX/MIN', 'qval')

out_coding$PV_location <- 'Coding'

absent <- c()
no_change <- c()
success <- c()
for(i in 1:nrow(PV_coding)){
  print('==========')
  out_coding$ID[i] <- PV_coding$PubMLST_id[i]
  
  if(!(PV_coding$PubMLST_id[i] %in% Exp$`1717.genes`)){
    absent <- c(absent, PV_coding$PubMLST_id[i])
    out_coding$Status[i] <- 'no_exp_data'
    next
  }
  
  out_coding$`MAX/MIN`[i] <- Exp[which(Exp$`1717.genes` == PV_coding$PubMLST_id[i]),14]
  out_coding$qval[i] <- Exp[which(Exp$`1717.genes` == PV_coding$PubMLST_id[i]),13]
  out_coding$Name[i] <- PV_coding$Name[i]
  out_coding$Function[i] <- PV_coding$Likely.Function[i]

  Exp_data <- Exp[which(Exp$`1717.genes` == PV_coding$PubMLST_id[i]),c(15:22)]
  colnames(Exp_data) <- c('27509', '27553', '28269', '28262', '53930', '28287', '53948', '53951')
  PV_data <- PV_coding[i, 17:24]
  colnames(PV_data) <- c('27509', '27553', '28262', '28269', '28287', '53930', '53948', '53951')
  combined_data <- rbind(PV_data, Exp_data)
  combined_data <- as.data.frame(t(combined_data))
  colnames(combined_data) <- c('state', 'expression')

  #Convert ON/OFF to a numeric 1 or 0
  for(j in 1:nrow(combined_data)){
    if(is.na(combined_data$state[j])){
      next
    }
    if(combined_data$state[j] == 'ON'){
      combined_data$state[j] <- 1
      next
    }
    if(combined_data$state[j] == 'OFF'){
      combined_data$state[j] <- 0
      next
    }
    if(combined_data$state[j] != 'ON' && combined_data$state[j] != 'OFF'){
      combined_data$state[j] <- NA
    }
  }

  #Skip genes where all coding loci are the same
  if(length(unique(na.omit(combined_data$state))) == 1){
    no_change <- c(no_change, PV_coding$PubMLST_id[i])
    out_coding$Status[i] <- 'no_state_change' 
    next
  }
  
  #Otherwise, run a wilcox test on combined_data, ensuring both columns are numeric
  combined_data$state <- as.numeric(combined_data$state)
  combined_data$expression <- as.numeric(combined_data$expression)
  out_wilcox_test$p.value <- NA
  try(out_wilcox_test <- t.test(combined_data$expression~combined_data$state))

  normality <- (shapiro.test(x=as.numeric(combined_data$expression)))
  
  if(normality$p.value > 0.05){
    print(paste('expression is normally distributed for', PV_coding$PubMLST_id[i]))
  } else{
    print(paste('expression is not normally distributed for', PV_coding$PubMLST_id[i]))
  }
  
  try(print(paste('t test p value is:',out_wilcox_test$p.value)))
  
  success <- c(success, PV_coding$PubMLST_id[i])
  out_coding$Status[i] <- 'state_change'
  try(out_coding$Pval[i] <- out_wilcox_test$p.value)
  #Plot ON/OFF against expression
  combined_data$state <- as.factor(combined_data$state)
  pdf(file = paste('~/Documents/PhD/PhD/PhasomeIt_data/Graphs/Coding_PV/', PV_coding$PubMLST_id[i], '.pdf', sep=''),
      width = 5, height = 6)
  p <- ggplot(na.omit(combined_data), aes(x=state, y=expression)) +
    geom_point() +
    ggtitle(paste(PV_coding$PubMLST_id[i])) +
    expand_limits(x=0, y =0) +
    coord_cartesian(expand = TRUE) +
    ggsignif::geom_signif(comparisons = list(c('0','1')), test = 't.test')
  print(p)
  dev.off()
}

absent #7 (6 unique) genes are missing from rnaseq
no_change #14 genes (13 unique) are in the same state in every isolate containing the tract
success #11 have a change in state

no_change <- c()
for(i in 1:nrow(PV_coding)){
  states <- unlist(PV_coding[i, 17:24])
  if(length(unique(na.omit(states))) == 1){
    no_change = c(no_change, PV_coding$PubMLST_id[i])
  }
}
length(unique(PV_IGR$PubMLST_id))

no_change <- c()
for(i in 1:nrow(PV_coding)){
  lengths <- unlist(PV_coding[i, 7:14])
  if(length(unique(na.omit(lengths))) == 1){
    no_change = c(no_change, PV_coding$PubMLST_id[i])
  }
}
length(unique(no_change))
length(unique(PV_IGR$PubMLST_id))

# 
# i=7
# PV_coding$PubMLST_id[i]
# Exp[which(Exp$`1717.genes` == 'NEIS0213'),15:22]
# 
# summary(lm(expression ~ state, data = x))
# combined_data[1,c(1:7)] <- 1
# combined_data[1,c(7:8)] <- 0
# x <- as.data.frame(t(combined_data))
# x$state <- as.numeric(x$state)
# x$expression <- as.numeric(x$expression)
# colnames(x) <- c('state', 'expression')
# 
# i=1
#For IGR tracts, perform linear regression
out_IGR <- as.data.frame(matrix(data = NA, nrow = nrow(PV_IGR), ncol = 9))
colnames(out_IGR) <- c('ID', 'Name', 'Function', 'PV_location', 'Status', 'Radj', 'Pval', 'MAX/MIN', 'qval')
#ID   IGR?    Status(absent,nochange,success)   Radj  Pval
out_IGR$PV_location <- 'IGR'

absent <- c()
no_change <- c()
success <- c()
for(i in 1:nrow(PV_IGR)){
  out_IGR$ID[i] <- PV_IGR$PubMLST_id[i]
  print(i)
  if(!(PV_IGR$PubMLST_id[i] %in% Exp$`1717.genes`) | is.na(PV_IGR$PubMLST_id[i])){
    absent <- c(absent, PV_IGR$PubMLST_id[i])
    out_IGR$Status[i] <- 'no_exp_data'
    next
  }
  
  out_IGR$`MAX/MIN`[i] <- Exp[which(Exp$`1717.genes` == PV_IGR$PubMLST_id[i]),14]
  out_IGR$qval[i] <- Exp[which(Exp$`1717.genes` == PV_IGR$PubMLST_id[i]),13]
  out_IGR$Name[i] <- PV_IGR$Name[i]
  out_IGR$Function[i] <- PV_IGR$Likely.Function[i]
  print('step 1')
  print(PV_IGR$PubMLST_id[i])
  
  Exp_data <- Exp[which(Exp$`1717.genes` == PV_IGR$PubMLST_id[i]),c(15:22)]
  colnames(Exp_data) <- c('27509', '27553', '28269', '28262', '53930', '28287', '53948', '53951')
  PV_data <- PV_IGR[i, 7:14]
  colnames(PV_data) <- c('27509', '27553', '28262', '28269', '28287', '53930', '53948', '53951')
  combined_data <- rbind(PV_data, Exp_data)
  combined_data <- as.data.frame(t(combined_data))
  colnames(combined_data) <- c('state', 'expression')
  
  print('step 2')
  print(combined_data)
  #Convert datapoints to numeric
  combined_data$state <- as.numeric(combined_data$state)
  combined_data$expression <- as.numeric(combined_data$expression)
  
  #If the tract is the same length in all isolates, skip
  if(length(unique(na.omit(combined_data$state))) == 1){
    no_change <- c(no_change, PV_IGR$PubMLST_id[i])
    out_IGR$Status[i] <- 'same_tract_length'
    next
  }
  
  model_out <- summary(lm(expression~state, data = na.omit(combined_data)))
  success <- c(success, PV_IGR$PubMLST_id[i])
  out_IGR$Status[i] <- 'change_tract_length'
  out_IGR$Radj[i] <- model_out$adj.r.squared
  out_IGR$Pval[i] <- model_out$coefficients[2,4]
  #Plot expression vs tract length
  pdf(file = paste('~/Documents/PhD/PhD/PhasomeIt_data/Graphs/IGR_PV/', PV_IGR$PubMLST_id[i], '.pdf', sep=''),
      width = 5, height = 6)
  p <- ggplot(na.omit(combined_data), aes(x=state, y=expression)) +
    geom_point() +
    ggtitle(paste(PV_IGR$PubMLST_id[i])) +
    expand_limits(y =0) +
    coord_cartesian(expand = TRUE)
  print(p)
  dev.off()
}

absent #3 entries missing from rnaseq data
no_change #6 entries with the same trct length
success #12 entries (11 unique) with a change


#from lm i need Radj and p value
#p value is summary(lm)$coefficiencts[2,4]

#Do i need to compare these loci to see
#get MAX/MIN logfold change and check the states of the gene in those two isolates, see if they differ. 
out_coding <- out_coding[order(out_coding$Pval),]
out_IGR <- out_IGR[order(out_IGR$Pval),]
#export out_coding and out_IGR to excel docs
#write.xlsx(out_coding,file = '~/Documents/PhD/PhD/PhasomeIt_data/coding_PV_stats.xlsx')
#write.xlsx(out_IGR, file = '~/Documents/PhD/PhD/PhasomeIt_data/IGR_PV_stats.xlsx')




#Create a summary table with the following headers:
#ID   Name    Function    Tract   Location    Repeat number range   n_isolates(contig breaks)   Exp(MAX/MIN)    P value

sumtable <- as.data.frame(matrix(data = NA, nrow = nrow(PVfilt), ncol = 9))
colnames(sumtable) <- c('ID', 'Name', 'Function', 'Tract', 'Location', 'Repeat number range', 'num isolates containing tract (contig breaks)', 'Exp (MAX/MIN)', 'p')
sumtable$ID <- PVfilt$PubMLST_id
sumtable$Name <- PVfilt$Name
sumtable$Function <- PVfilt$Likely.Function
sumtable$Tract <- PVfilt$Tract

for(i in 1:nrow(PVfilt)){
  print(i)
  #Obtain location of the tract
  Loc <- PVfilt[i, 16:24]
  Loc <- unique(na.omit(unlist(Loc)))
  if(Loc[1] == 'IGR'){
    sumtable$Location[i] <- 'IGR'
  }
  if(Loc[1] != 'IGR'){
    sumtable$Location[i] <- 'CDS'
  }
  
  #find min and max repeat length
  repeats <- as.vector(na.omit(unlist(PVfilt[i, 7:14])))
  entry <- paste(min(repeats),max(repeats), sep =' - ')
  if(min(repeats) == max(repeats)){
    entry <- min(repeats)
  }
  sumtable$`Repeat number range`[i] <- entry
  
  #Find number of isolates that contain the complete tract
  #how many were indeterminate due to contig breaks needs manual checking
  tracts <- na.omit(unlist(PVfilt[i, 17:24]))
  sumtable$`num isolates containing tract (contig breaks)`[i] <- length(tracts)
  
  #Find the expression for the max and min isolates at that locus
  if(PVfilt$PubMLST_id[i] %in% Exp$`1717.genes`){
    
    sumtable$`Exp (MAX/MIN)`[i] <- Exp[which(Exp$`1717.genes` == PVfilt$PubMLST_id[i]), 'log2(EXP/Min)']
  }
  
  #get the p value for the earlier signfiicance testing 
  if(sumtable$Location[i] == 'IGR'){
    sumtable$p[i] <- out_IGR[which(out_IGR$ID == sumtable$ID[i]),'Pval']
  }
  if(sumtable$Location[i] == 'CDS'){
    sumtable$p[i] <- out_coding[which(out_coding$ID == sumtable$ID[i]),'Pval']
  }
  
}

sumtable <- sumtable[order(sumtable$Location),]
write.xlsx(x = sumtable, file = '~/Documents/PhD/PhD/PhasomeIt_data/summary_tabletemp.xlsx')


##Check which PV genes are significant hits in rnaseq
rna_degs <- c()
for(i in PV$PubMLST_id){
  
}
PV_rna_degs <- Exp[which(Exp$`1717.genes` %in% PV$PubMLST_id), c(9,13,14)]
missing <- PV$PubMLST_id[which(!(PV$PubMLST_id %in% PV_rna_degs$`1717.genes`))]
missing
PV_rna_degs$X13 <- as.numeric(PV_rna_degs$X13)
PV_rna_degs <- PV_rna_degs[which(PV_rna_degs$X13 < 0.05),]
PV_rna_degs <- PV_rna_degs[which(PV_rna_degs$`log2(EXP/Min)` >= 1),]

degs <- PV[which(PV$PubMLST_id %in% PV_rna_degs$`1717.genes`),]



###For each PV igr gene check if significantly different isolates have different igr repeats 
#Remake rna dataframe and sort to have just our columns of interest
rna <- read.xlsx("~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx", sheet = 8)

#Only columns needed are the gene names in NMB and NEIS, their products, expression values, and pairwise qvalues
rna <- rna[, c(8:10,12, 15:22, 101:156) ]

colnames(rna) <- c("CDS_code", "NEIS_id", "NMB_id", "Product", 
                   "exp_20578", "exp_22689", "exp_N188.1", "exp_N222.1", "exp_N222.2", "exp_N445.1", "exp_N459.3", "exp_N459.6",
                   "q20578_22689", "q20578_N188.1", "q20578_N222.1", "q20578_N222.2", "q20578_N445.1", "q20578_N459.3", "q20578_N459.6",
                   "q22689_20578", "q22689_N188.1", "q22689_N222.1", "q22689_N222.2", "q22689_N445.1", "q22689_N459.3", "q22689_N459.6",
                   "qN188.1_20578", "qN188.1_22689", "qN188.1_N222.1", "qN188.1_N222.2", "qN188.1_N445.1", "qN188.1_N459.3", "qN188.1_N459.6",
                   "qN222.1_20578", "qN222.1_22689", "qN222.1_N188.1", "qN222.1_N222.2", "qN222.1_N445.1", "qN222.1_N459.3", "qN222.1_N459.6",
                   "qN222.2_20578", "qN222.2_22689", "qN222.2_N188.1", "qN222.2_N222.1", "qN222.2_N445.1", "qN222.2_N459.3", "qN222.2_N459.6",
                   "qN445.1_20578", "qN445.1_22689", "qN445.1_N188.1", "qN445.1_N222.1", "qN445.1_N222.2", "qN445.1_N459.3", "qN445.1_N459.6",
                   "qN459.3_20578", "qN459.3_22689", "qN459.3_N188.1", "qN459.3_N222.1", "qN459.3_N222.2", "qN459.3_N445.1", "qN459.3_N459.6",
                   "qN459.6_20578", "qN459.6_22689", "qN459.6_N188.1", "qN459.6_N222.1", "qN459.6_N222.2", "qN459.6_N445.1", "qN459.6_N459.3"
)
#Remove first row
rna = rna[-c(1),]
rna[,5:68] <- as.numeric(rna[,5:68])

cols <- colnames(rna)[5:68]
rna[cols] <- sapply(rna[cols], as.numeric)


check_repeats <- function(locus){
  
}

isolates <- c('i1', 'i2')
for(i in 1:length(isolates)){
  print(i)
}
