operondat <- read.xlsx('~/Documents/PhD/PhD/operon_mapper_res/opdata.xlsx')
rna <- read.xlsx('~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-04092024.xlsx', sheet = 8)

operondat <- operondat[which(!is.na(operondat$pubmlst_id)),]
#remove whitespace in rna
rna$`1717.genes`<-str_trim(rna$`1717.genes`)
rna$X13 <- as.numeric(rna$X13)
rna$`log2(EXP/Min)` <- as.numeric(rna$`log2(EXP/Min)`)

#how many genes in operons?, excluding na entries
nrow(operondat) #1120

#how many genes in operons are significant?
rna <- rna[which(rna$`1717.genes` %in% operondat$pubmlst_id & !is.na(rna$`1717.genes`)),]
nrow(rna) #945 after removing na pubmlst entries

#how many of those are significant
sig <- rna[which(as.numeric(rna$X13) < 0.05 & rna$`log2(EXP/Min)` >= 1),]
nrow(sig) #239

239/945 #25.3%

#How many genes not in operons, excluding na entries

operondat <- read.xlsx('~/Documents/PhD/PhD/operon_mapper_res/opdata.xlsx')
rna <- read.xlsx('~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-04092024.xlsx', sheet = 8)

operondat <- operondat[which(!is.na(operondat$pubmlst_id)),]
#remove whitespace in rna
rna$`1717.genes`<-str_trim(rna$`1717.genes`)
rna$X13 <- as.numeric(rna$X13)
rna$`log2(EXP/Min)` <- as.numeric(rna$`log2(EXP/Min)`)

nonoperondat <- rna[which(!(rna$`1717.genes` %in% operondat$pubmlst_id)),]
nrow(nonoperondat)
nonoperondat <- nonoperondat[which(!is.na(nonoperondat$`1717.genes`)),]
nrow(nonoperondat) #779 non operon genes excluding na

#how many nonoperon non NA genes are significant?
sig <- nonoperondat[which(as.numeric(nonoperondat$X13) < 0.05 & nonoperondat$`log2(EXP/Min)` >= 1),]
nrow(sig) #184
184/779 #23.6%
