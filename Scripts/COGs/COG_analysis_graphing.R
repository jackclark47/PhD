#Script to get COG tables for each isolate, summarising EGGNOG outputs
library(openxlsx)
library(ggplot2)
library(gplots)
library(corrplot)

setwd("~/Documents/PhD/RNA_IGR/eggNOG_annotations/Isolates")

COGsum <- function(isolate, append){
  data <- read.xlsx(paste(isolate, '.annotations.xlsx', sep = ''), startRow = 3)[,7]
  annotated <- length(data)
  return(table(data))
}

i1 <- COGsum('27509', TRUE)
i2 <- COGsum('27553', TRUE)
i3 <- COGsum('28262', TRUE)
i4 <- COGsum('28269', TRUE)
i5 <- COGsum('28287', TRUE)
i6 <- COGsum('53930', TRUE)
i7 <- COGsum('53948', TRUE)
i8 <- COGsum('53951', TRUE)

df <- as.data.frame(cbind(i1,i2,i3,i4,i5,i6,i7,i8))
colnames(df) <- c('27509', '27553', '28262', '28287', '53930', '53948', '53951')
setwd("~/Documents/PhD/RNA_IGR/eggNOG_annotations/")
# write.xlsx(df, "eggNOG_COG_summaries_new.xlsx")


# write.xlsx(table(read.xlsx('~/Downloads/eggNOG_annotations/all_coding_loci.annotations.xlsx', startRow = 3)[,7]), 'all_loci_COG_summaries.xlsx')
# write.xlsx(table(read.xlsx('~/Downloads/eggNOG_annotations/N222.2.annotations.xlsx', startRow = 3)[,7]), 'N222.2.annotations.xlsx')



#Section to compare the all_loci annotations to those in the rna expression dataset, and laater to the igr dataset
#load rna loci
rna <- read.xlsx("~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-04092024.xlsx", sheet=8)
rnanonsig <- rna[(which(as.numeric(rna$X13) < 0.05 | rna$`log2(EXP/Min)` < 1)), ]
rna_underlog2 <- rna[which(rna$`log2(EXP/Min)` < 1),]
rna <- rna[which(rna$`log2(EXP/Min)` >=1),]
rnalog4 <- rna[which(rna$`log2(EXP/Min)` >=2),]
rnasig <- rna[which(as.numeric(rna$X13) < 0.05),]

#load distributions of cogs in all loci
cogs <- read.xlsx('~/Documents/PhD/RNA_IGR/eggNOG_annotations/all_coding_loci.annotations.xlsx', startRow = 3)[1:1765,]

#match loci in the two lists using log2fold change >1 and >2 as cutoffs

cogslog2 <- cogs[which(cogs$query %in% rna$`1717.genes`),] 
length(unique(rna$`1717.genes`)) #780 loci, of which 546 have cog information
cogslog4 <- cogs[which(cogs$query %in% rnalog4$`1717.genes`),] 
length(unique(rnalog4$`1717.genes`)) #165 loci, of which 115 have cog information
cogsunderlog2 <- cogs[which(cogs$query %in% rna_underlog2$`1717.genes`),] 

cogssig <- cogs[which(cogs$query %in% rnasig$`1717.genes`),] 
cogsnonsig <- cogs[which(cogs$query %in% rnanonsig$`1717.genes`),] 

#get table of COG categories for those two cutoffs
dataset <- list('log2>1' = table(cogslog2$COG_category), 'log2>2' = table(cogslog4$COG_category), 'under_log2' = table(cogsunderlog2$COG_category))
#write.xlsx(dataset, file = 'rna_cogstest.xlsx')

dataset <- list('sig' = table(cogssig$COG_category), 'nonsig' = table(cogsnonsig$COG_category))
#write.xlsx(dataset, file = 'sig_rna_cogs.xlsx')



#now for igr categories
#load lists of loci present in each category
upvar <- read.xlsx("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_all.xlsx", sheet = 7)
downvar <- read.xlsx("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_all.xlsx", sheet = 4)
highvar <- read.xlsx("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_all.xlsx", sheet = 5)
lowvar <- read.xlsx("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_all.xlsx", sheet = 6)
novar <- read.xlsx("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_all.xlsx", sheet = 2)
var <- read.xlsx("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/Correct/igrs_all.xlsx", sheet = 3)
#what are the cog distributions in each category?
upvarcogs <- cogs[which(cogs$query %in% upvar$Locus),][,c(1,7)]
downvarcogs <- cogs[which(cogs$query %in% downvar$Locus),][,c(1,7)]
highvarcogs <- cogs[which(cogs$query %in% highvar$Locus),][,c(1,7)]
lowvarcogs <- cogs[which(cogs$query %in% lowvar$Locus),][,c(1,7)]
novarcogs <- cogs[which(cogs$query %in% novar$Locus),][,c(1,7)]
varcogs <- cogs[which(cogs$query %in% var$Locus),][,c(1,7)]

dataset <- list('upvar' = table(upvarcogs$COG_category), 'downvar' = table(downvarcogs$COG_category) , 'highvar' = table(highvarcogs$COG_category), 'lowvar' =table(lowvarcogs$COG_category), 'novar' = table(novarcogs$COG_category), 'var' = table(varcogs$COG_category))
#write.xlsx(dataset, file = 'igr_cogs.xlsx')

#get balloon plots of cog categories by expression and by pi, using chi squared to test for significance as we have counts data


#first need to make contingency tables for each dataset of the format:
#       All_loci    log2    log4
# A         100       78      77
# B         57
# C
# D
#etc etc

categories <- c("Unassigned", "C: Energy production/conversion", "D: Cell cycle, division", "E: AA transport/metabolism", "F: Nucleotide transport/metabolism",
                "G: Carbohydrate transport/metabolism", "H: Coenzyme transport/metabolism", "I: Lipid transport/metabolism", "J: Translation, ribosomes, biogenesis",
                "K: Transcription", "L: Replication, Recombination, Repair", "M: Cell wall/membrane/envelope biogenesis", "N: Cell motility", "O: PTM, protein turnover, chaperones",
                "P: Inorganic ion transport/metabolism", "Q: Secondary metabolite transport/metabolism/synthesis", "S: Unknown", "T: Signal transduction mechanisms",
                "U: Intracellular trafficking, secretion, vesicular transport", "V: Defense mechanisms")

#for rna data
rna_counts <- read.xlsx("~/Documents/PhD/RNA_IGR/eggNOG_annotations/rna_cogstest.xlsx")[1:22,]
colnames(rna_counts) <- c('COG_category', '2fold', '4fold', 'under_2fold')
rownames(rna_counts) <- rna_counts$COG_category
rna_counts <- rna_counts[,-1]
#for igr data
igr_counts <- read.xlsx("~/Documents/PhD/RNA_IGR/eggNOG_annotations/igr_cogs.xlsx")[1:22,c(3:9)]
colnames(igr_counts) <- c('COG_category', 'upvar', 'downvar', 'highvar', 'lowvar', 'novar', 'var')
rownames(igr_counts) <- igr_counts$COG_category
igr_counts <- igr_counts[,-1]

sig_counts <- read.xlsx("~/Documents/PhD/RNA_IGR/eggNOG_annotations/sig_rna_cogs.xlsx")[1:20,c(1:3)]
colnames(sig_counts) <- c('COG_category', 'sig', 'nonsig')
rownames(sig_counts) <- sig_counts$COG_category
sig_counts <- sig_counts[,-1]
#need to exclude data where the expected is less than 5 counts
igr_counts <- igr_counts[-c(2,22),]
rna_counts <- rna_counts[-c(2,22),]
categories <- c("Unassigned", "C: Energy production/conversion", "D: Cell cycle, division", "E: AA transport/metabolism", "F: Nucleotide transport/metabolism",
                "G: Carbohydrate transport/metabolism", "H: Coenzyme transport/metabolism", "I: Lipid transport/metabolism", "J: Translation, ribosomes, biogenesis",
                "K: Transcription", "L: Replication, Recombination, Repair", "M: Cell wall/membrane/envelope biogenesis", "N: Cell motility", "O: PTM, protein turnover, chaperones",
                "P: Inorganic ion transport/metabolism", "Q: Secondary metabolite transport/metabolism/synthesis", "S: Unknown", "T: Signal transduction mechanisms",
                "U: Intracellular trafficking, secretion, vesicular transport", "V: Defense mechanisms")
rownames(igr_counts) <- categories
rownames(rna_counts) <- categories
rownames(sig_counts) <- categories

#compute chisquared tests
rna_chisq <- chisq.test(rna_counts)
rna_chisq$p.value #significant between 2 fold and below 2 fold. barely not significant between all 3
corrplot(rna_chisq$residuals, method = 'color', is.cor = FALSE, tl.cex = 0.9, tl.col = 'black', 
         cl.ratio = 0.5,
         cl.align.text = 'l',
         cl.offset = 0.2)

igr_chisq <- chisq.test(igr_counts[,c(1:5)])
igr_chisq$p.value
corrplot(igr_chisq$residuals, method = 'color',is.cor = FALSE, tl.cex = 1, tl.col = 'black',
         cl.ratio = 0.3,
         cl.align.text = 'l',
         cl.offset = 0.2)


var_chisq <- chisq.test(igr_counts[,c(5:6)])
var_chisq$p.value
corrplot(var_chisq$residuals, method = 'color',is.cor = FALSE, tl.cex = 1, tl.col = 'black',
         cl.ratio = 0.8,
         cl.align.text = 'l',
         cl.offset = 0.2)


sig_chisq <- chisq.test(sig_counts)
sig_chisq$p.value
corrplot(sig_chisq$residuals, method = 'color',is.cor = FALSE, tl.cex = 1, tl.col = 'black',
         cl.ratio = 0.8,
         cl.align.text = 'l',
         cl.offset = 0.2)


rna_chisq$observed
round(rna_chisq$expected,2)

contrib_rna <- 100*rna_chisq$residuals^2/rna_chisq$statistic
contrib_rna <- round(contrib_rna,3)
corrplot(contrib_rna, is.corr = F)


contrib_igr <- 100*igr_chisq$residuals^2/igr_chisq$statistic
contrib_igr <- round(contrib_igr,3)
corrplot(contrib_igr, is.corr = F)






#Repeat for each isolate to see if there are differences between them

#load isolate spreadsheet
isolate_cogs <- read.xlsx("~/Documents/PhD/RNA_IGR/eggNOG_annotations/Old/eggNOG_COG_summaries.xlsx")[1:22,9:16]
isolate_cogs <- isolate_cogs[-c(2,22),]
rownames(isolate_cogs) <- categories
isolate_cogs <- isolate_cogs[,-c(1)]
isolate_cogs$`27509` <- as.numeric(isolate_cogs$`27509`)
isolate_cogs$`27553` <- as.numeric(isolate_cogs$`27553`)
isolate_cogs$`28262` <- as.numeric(isolate_cogs$`28262`)
isolate_cogs$`28287` <- as.numeric(isolate_cogs$`28287`)
isolate_cogs$`53930` <- as.numeric(isolate_cogs$`53930`)
isolate_cogs$`53948` <- as.numeric(isolate_cogs$`53948`)
isolate_cogs$`53951` <- as.numeric(isolate_cogs$`53951`)

isolate_chisq <- chisq.test(isolate_cogs)
isolate_chisq$p.value
corrplot(isolate_chisq$residuals, is.cor = FALSE, tl.cex = 1, tl.col = 'black')

contrib_isolate <- 100*isolate_chisq$residuals^2/isolate_chisq$statistic
contrib_isolate <- round(contrib_isolate,3)
corrplot(contrib_isolate, is.corr = F)
