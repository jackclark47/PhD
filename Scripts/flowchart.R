#Read in the datasets
Exp <- read.xlsx('~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-16082022.xlsx', sheet=8)
novar <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 2)
var <- read.xlsx('~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx', sheet = 3)
PV <- read.xlsx('~/Documents/PhD/PhD/PhasomeIt_data/phasomeit_out_RNAcheck.xlsx', sheet = 5)
opdata <- read.xlsx('~/Documents/PhD/PhD/operon_mapper_res/opdata.xlsx', sheet = 6)
op_summary <- read.xlsx('~/Documents/PhD/PhD/operon_mapper_res/operon_summary.xlsx')
#Preprocess
PV <- PV[which(!is.na(PV$PubMLST_id)),]
Exp$`1717.genes`<-str_trim(Exp$`1717.genes`)
Exp$X13 <- as.numeric(Exp$X13)
Exp$`log2(EXP/Min)` <- as.numeric(Exp$`log2(EXP/Min)`)


#Get sig genes in Exp and those with no annotation
degs <- Exp[which(Exp$`log2(EXP/Min)` >= 1),]
degs <- degs[which(degs$X13 < 0.05),]
missing <- degs[which(is.na(degs$`1717.genes`)),]


#Get sig var and no var. Identify number of significant variable genes unannotated in exp data
sig_var <- degs[which(degs$`1717.genes` %in% var$Locus),]
sig_novar <- degs[which(degs$`1717.genes` %in% novar$Locus),]
missing_sigvar <- degs[which(!(degs$`1717.genes` %in% sig_var$`1717.genes`)),]
missing_sigvar <- degs[which(!(missing_sigvar$`1717.genes` %in% sig_novar$`1717.genes`)),]

nrow(sig_var) #177
nrow(sig_novar) #122
nrow(missing_sigvar) #90


#Get nonsig var and no var. Identify number of nonsignificant variable genes unannotated in exp data
ns_genes <- Exp[which(!(Exp$`1717.genes` %in% degs$`1717.genes`)),]
nsig_var <- ns_genes[which(ns_genes$`1717.genes` %in% var$Locus),]
nsig_novar <- ns_genes[which(ns_genes$`1717.genes` %in% novar$Locus),]
nsmissing_var <- ns_genes[which(!(ns_genes$`1717.genes` %in% nsig_var$`1717.genes`)),]
nsmissing_var <- nsmissing_var[which(!(nsmissing_var$`1717.genes` %in% nsig_novar$`1717.genes`)),]

nrow(nsig_var) #681
nrow(nsig_novar) #725
nrow(nsmissing_var) #13

681 + 725 + 13 == nrow(ns_genes)


#Get number of significant PV genes without igrs. Am including operon based PV genes here
sig_PV_novar <- sig_novar[which(sig_novar$`1717.genes` %in% PV$PubMLST_id),]
nrow(sig_PV_novar) #just 1 
#How many sig PV genes with igrs are there?
sig_PV_var <- sig_var[which(sig_var$`1717.genes` %in% PV$PubMLST_id),]
nrow(sig_PV_var) #8 

sig_PV <- degs[which(degs$`1717.genes` %in% PV$PubMLST_id),]
nrow(sig_PV)
sig_PV[which(sig_PV$`1717.genes` %in% sig_var$`1717.genes`),]


#How many IGRPV genes are significant in the rna?
igr_pv_rows <- c()
cds_pv_rows <- c()
for(i in 1:nrow(PV)){
  states <- as.vector(na.omit(unique(unlist(as.vector(PV[i,31:38])))))
  print(states)[1]
  try(
  if(states[1] == 'IGR'){
    igr_pv_rows <- c(igr_pv_rows,i)
  })
  try(
    if(states[1] == 'ON' | states[1] == 'OFF'){
      cds_pv_rows <- c(cds_pv_rows, i)
    }
  )
}

#How many CDSPV genes are significnat in the rna
igr_pvs <- PV[igr_pv_rows,]
cds_pvs <- PV[cds_pv_rows,]
sig_igr_pvs <- igr_pvs[which(igr_pvs$PubMLST_id %in% degs$`1717.genes`),]
sig_cds_pvs <- cds_pvs[which(cds_pvs$PubMLST_id %in% degs$`1717.genes`),]
nrow(sig_igr_pvs) #3
nrow(sig_cds_pvs) #7


#Now get number of significant genes without igrs in operons
operonic <- opdata[which(opdata$pubmlst_id %in% sig_novar$`1717.genes`),]
nrow(operonic)
operonic <- operonic[which(operonic$sig_gene == 1),]
nrow(operonic)

#How many sig genes are in operons with igr variation?
variable_opsum <- op_summary[which(op_summary$variable==1),]
variable_opdata <- opdata[which(opdata$Operon %in% variable_opsum$Operon),]
sig_igr_operonic <- operonic[which(operonic$pubmlst_id %in% variable_opdata$pubmlst_id),]
nrow(sig_igr_operonic)

nrow(op_summary[which(op_summary$variable==1),])
nrow(op_summary[which(op_summary$variable==0),])

nrow(opdata)
operonic_degs <- opdata[which(opdata$pubmlst_id %in% degs$`1717.genes`),]
nrow(operonic_degs)
nrow(operonic_degs[which(operonic_degs$pubmlst_id %in% novar$Locus),])




#Get info for VarGroups


