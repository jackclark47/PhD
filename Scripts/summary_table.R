#This script collates data from all the data sets - PV, operon, RNA, IGR, CDS, and produces a summary table 
#The table contains info on which groups each locus belongs to, to make it easier to quickly get useful summaries

#Load in datasets
Exp <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx', sheet=8)[-1,]
novar <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 2)
novar <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 2)
upvar <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 6)
downvar <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 5)
highvar <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 7)
lowvar <- read.xlsx('~/Documents/PhD/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 4)
PV <- read.xlsx('~/Documents/PhD/PhD/PhasomeIt_data/phasomeit_out_RNAcheck.xlsx', sheet = 5)
opdata <- read.xlsx('~/Documents/PhD/PhD/operon_mapper_res/opdata.xlsx', sheet = 6)
op_summary <- read.xlsx('~/Documents/PhD/PhD/operon_mapper_res/operon_summary.xlsx')

#Still need:
#CDS
#Regulon
#Phasome

#Preprocess
PV <- PV[which(!is.na(PV$PubMLST_id)),]
Exp$`1717.genes`<-str_trim(Exp$`1717.genes`)
Exp$X13 <- as.numeric(Exp$X13)
Exp$`log2(EXP/Min)` <- as.numeric(Exp$`log2(EXP/Min)`)


#Initialise the dataframe and add locus names from rna data
summary <- as.data.frame(matrix(data = NA, nrow = nrow(Exp), ncol = 9))
colnames(summary) <- c('pubmlst_id', 'rna_sig', 'var_igr', 'pv', 'pv_change', 'operon', 'var_cds', 'regulon', 'phasome')
summary$pubmlst_id <- Exp$`1717.genes`

#get significant genes and assign them to the ids
summary$rna_sig <- 0
summary$rna_sig[which(Exp$X13 < 0.05 & Exp$`log2(EXP/Min)` >=1)] <- 1

#Add PV data - pv column is just saying if the gene is pv or not. pv_chaneg is if theres a change in state (CDS) or repeat tract length (igr)
summary$pv <- 0
summary$pv_change <- 0
summary$pv[which(summary$pubmlst_id %in% PV$PubMLST_id)] <- 1
#summary$pv_change[which(summary$pubmlst_id %in% PV$PubMLST_id & (PV$`Change.in.state.or.frame.shift?` == 1 | PV$`Change.in.state.or.frame.shift?` == '?'))] <- 1
ids <- PV$PubMLST_id[which(PV$`Change.in.state.or.frame.shift?` == 1 | PV$`Change.in.state.or.frame.shift?` == '?')]
summary$pv_change[which(summary$pubmlst_id %in% ids)] <- 1

#Add igr data
summary$var_igr[which(summary$pubmlst_id %in% upvar$Locus)] <- 'up_var'
summary$var_igr[which(summary$pubmlst_id %in% novar$Locus)] <- 'no_var'
summary$var_igr[which(summary$pubmlst_id %in% downvar$Locus)] <- 'down_var'
summary$var_igr[which(summary$pubmlst_id %in% highvar$Locus)] <- 'high_var'
summary$var_igr[which(summary$pubmlst_id %in% lowvar$Locus)] <- 'low_var'

#add operon data
sig_ops <- opdata[which(opdata$sig_operon == 1),]
summary$operon <- 0
summary$operon[which(summary$pubmlst_id %in% opdata$pubmlst_id)] <- 'operonic'

summary$operon[which(summary$pubmlst_id %in% sig_ops$pubmlst_id)] <- 'sig_operonic'
summary$operon[which(is.na(summary$pubmlst_id))] <- 0

write.xlsx(summary, '~/Documents/PhD/PhD/all_summary.xlsx')
#Add CDS data


#Add phasome 



#Add regulon

