#This is a script to analyse RNA expression data and IGR variation in RNA hits.
#First I obtain a list of loci which display a >log2 fold change in expression
#Then using previous data on IGR variation in those loci, I group them based on no IGR variation, or IGR variation
#The IGR variaiton group is further split into whether it is the up or downstream IGR, and the degree of variation
#I then run stats across those groupings to see if there are more expression changes in the IGR variants than the no variation dataset
#Then I visualise.

#Load modules
library(ggplot2)
library(openxlsx)
library(ggpubr)
library(dplyr)
library(stringr)
setwd("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data")

#Read in data
log2hits <- read.xlsx("~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-04092024.xlsx", sheet = 8) #Change this depending on if you want log(EXP/MIN > 2 (sheet 10), include reads < 50 (sheet 9), or all loci (sheet 8))
noVar <- read.xlsx("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx", sheet = 2)
upVar <- read.xlsx("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx", sheet = 6)
downVar <- read.xlsx("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx", sheet = 5)
lowVar <- read.xlsx("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx", sheet = 4)
highVar <- read.xlsx("~/Documents/PhD/RNA_IGR/Isolate_Igr_Data/igr_new_11012024.xlsx", sheet = 7)

#noVar, loci with no IGR variation
#upVar, loci with variation only in their upstream IGR
#downVar, the inverse of upVar
#lowVar, loci with variation in both up and downstream IGRs, but a score lower than 24 when summing allele rows
#highVar, loci with variaiton in both IGRs and a score of 24 or greater
#log2hits are loci from RNAseq which hace more than 50 maxreads and have been filtered to exclude things like predicted RNAs, ribosomal proteins


#Optional filtering to obtain only those with q values < 0.05 i.e. the significant hits
log2hits$X13 <- as.numeric(log2hits$X13)
log2hits <- log2hits[which(log2hits$X13 < 0.05),]
log2hits <- log2hits[which(log2hits$`log2(EXP/Min)` >= 1),]




#Remove whitespaces in locus names in the rna dataset
for(i in 1:length(log2hits$`1717.genes`)){
  log2hits$`1717.genes`[i] <- str_trim(log2hits$`1717.genes`[i])
}





###Obtain list of significant hits from RNAseq
#Reduce number of columns and remove first row
log2hitsfilt <- log2hits[-1,-c(1:8, 10, 11, 36:107)]

#Remove entries which don't have an NEIS code
log2hitsfilt <- log2hitsfilt[!is.na(log2hitsfilt$`1717.genes`),] #548 have no NEIS code

#Correct column names 
colnames(log2hitsfilt) <- c("PubMLST_id", "Product", "Min_qvalue", "Max_log2", "exp_20578_1", "exp_22689_1", "exp_N188_1",
                            "exp_N222_1", "exp_N222_2", "exp_N445_1", "exp_N459_3", "exp_N459_6", "Min_exp", "Min_ID",
                            "log2_20578_1", "log2_22689_1", "log2_N188_1", "log2_N222_1", "log2_N222_2", "log2_N445_1", "log2_N459_3",
                            "log2_N459_6", "Max_log2", "Min_reads", "Max_exp")

#Create empty data frame of the data thatll be used for plotting
plotdata <- as.data.frame(matrix(nrow=1723,ncol=4)) #524 or 1723
alleleplotdata <- as.data.frame(matrix(nrow=1723,ncol=4))
colnames(plotdata) <- c('Locus', 'log2fold.change', 'VarGroup', 'Var')
colnames(alleleplotdata) <- c('Locus', 'log2fold.change', 'VarGroup', 'Allele.count')

###Loop through list log2hitsfilt and sort into groups
#Initialise lists
nomatch <- c()
noVar_hits <- c()
upVar_hits <- c()
downVar_hits <- c()
lowVar_hits <- c()
highVar_hits <- c()
Var_hits <- c() #Store every hit which has variation in an igr somewhere
multihit <- c() #Store loci appearing in multiple of the above lists, except Var_hits and nomatch

for(i in 1:nrow(log2hitsfilt)){
  count = 0 #Count number of times that locus is added to a list
  id <- log2hitsfilt$PubMLST_id[i]
  plotdata$Locus[i] <- id
  plotdata$log2fold.change[i] <- log2hitsfilt$Max_log2[i]
  alleleplotdata$Locus[i] <- id
  alleleplotdata$log2fold.change[i] <- log2hitsfilt$Max_log2[i]
  #Check no var
  if(id %in% noVar$Locus){
    noVar_hits <- c(noVar_hits, id)
    count = count + 1
    plotdata$VarGroup[i] <- "NoVar"
    plotdata$Var[i] <- "NoVar"
    alleleplotdata$VarGroup[i] <- "NoVar"
    alleleplotdata$Allele.count[i] <- 1
  }
  #Check upVar
  if(id %in% upVar$Locus){
    upVar_hits <- c(upVar_hits, id)
    Var_hits <- c(Var_hits, id)
    count = count + 1
    plotdata$VarGroup[i] <- "upVar"
    plotdata$Var[i] <- "Var"
    alleleplotdata$VarGroup[i] <- "upVar"
    alleleplotdata$Allele.count[i] <- upVar[upVar$Locus==id,]$Combined
  }
  #Check downVar
  if(id %in% downVar$Locus){
    downVar_hits <- c(downVar_hits, id)
    Var_hits <- c(Var_hits, id)
    count = count + 1
    plotdata$VarGroup[i] <- "downVar"
    plotdata$Var[i] <- "Var"
    alleleplotdata$VarGroup[i] <- "downVar"
    alleleplotdata$Allele.count[i] <- downVar[downVar$Locus==id,]$Combined
  }
  #Check lowVar
  if(id %in% lowVar$Locus){
    lowVar_hits <- c(lowVar_hits, id)
    Var_hits <- c(Var_hits, id)
    count = count + 1 
    plotdata$VarGroup[i] <- "lowVar"
    plotdata$Var[i] <- "Var"
    alleleplotdata$VarGroup[i] <- "lowVar"
    alleleplotdata$Allele.count[i] <- lowVar[lowVar$Locus==id,]$Combined
  }
  #Check highVar
  if(id %in% highVar$Locus){
    highVar_hits <- c(highVar_hits, id)
    Var_hits <- c(Var_hits, id)
    count = count + 1
    plotdata$VarGroup[i] <- "highVar"
    plotdata$Var[i] <- "Var"
    alleleplotdata$VarGroup[i] <- "highVar"
    alleleplotdata$Allele.count[i] <- highVar[highVar$Locus==id,]$Combined
  }
  #if in none add to list of anomalies
  if(count == 0){
    nomatch <- c(nomatch, id)
    plotdata$VarGroup[i] <- "NoMatch"
    alleleplotdata$VarGroup[i] <- "NoMatch"
  }
  if(count > 1){
    multihit <- c(multihit, id)
  }
}

table(alleleplotdata$VarGroup)
highVar$Locus
highVar[highVar$Locus==id,]$Combined
id='NEIS2198'

#Check the core groups are the same size as the original dataset
length(Var_hits) + length(noVar_hits) + length(nomatch) == length(log2hitsfilt$PubMLST_id)
#TRUE
length(multihit)
length(nomatch) #120 significant rna hits had no igr variation or were not contained in the dataset



###Finding basic stats
#How many loci were looked at in total in RNAseq? 
nloci <- 2271 #423 or 2271
#How many loci are variable in at least one IGR?
#How many aren't variable at all?
RNA_loci <- read.xlsx("~/Documents/PhD/RNA_IGR/11-7_run_5_transcripts-cdb-04092024.xlsx", sheet = 8)[-c(1),]
for(i in 1:length(RNA_loci$`1717.genes`)){
  RNA_loci$`1717.genes`[i] <- str_trim(RNA_loci$`1717.genes`[i])
}

#RNA_loci <- log2hits #enable if using only significant rnaseq loci
Varcount = 0 #Initialise counts of each type of locus for the full dataset, not just significant hits
noVarcount = 0
NAcount = 0 
lowVarcount = 0
highVarcount = 0
upVarcount = 0
downVarcount = 0
#Loop that adds to the counts depending on where that rnaseq locus appears in the igr dataset
for(id in 1:nrow(RNA_loci)){
  if(is.na(RNA_loci$`1717.genes`[id])){
    NAcount = NAcount + 1 
    next
  }
  if(RNA_loci$`1717.genes`[id] %in% highVar$Locus | RNA_loci$`1717.genes`[id] %in% lowVar$Locus | RNA_loci$`1717.genes`[id] %in% upVar$Locus | RNA_loci$`1717.genes`[id] %in% downVar$Locus){
    Varcount = Varcount +1
    if(RNA_loci$`1717.genes`[id] %in% highVar$Locus){
      highVarcount = highVarcount + 1
    }
    if(RNA_loci$`1717.genes`[id] %in% lowVar$Locus){
      lowVarcount = lowVarcount + 1
    }
    if(RNA_loci$`1717.genes`[id] %in% upVar$Locus){
      upVarcount = upVarcount + 1
    }
    if(RNA_loci$`1717.genes`[id] %in% downVar$Locus){
      downVarcount = downVarcount + 1
    }
    next
  }
  else{
    noVarcount = noVarcount + 1
  }
}      #old dataset.      #new dataset
NAcount #554              554
Varcount #863             719
noVarcount #854           998
lowVarcount #345          290
highVarcount #99          43
upVarcount # 218          194
downVarcount #201         192
NAcount + noVarcount + Varcount == nrow(RNA_loci)
lowVarcount + highVarcount + upVarcount + downVarcount == Varcount
Varcounts_all <- as.data.frame(list(NAcount, Varcount, noVarcount, lowVarcount, highVarcount, upVarcount, downVarcount),
                           col.names = c("NA", "Var", "noVar", "lowVar", "highVar", "upVar", "downVar"))
noVarcount + Varcount == nrow(log2hitsfilt)
###Now to make a plot of boxplots per number of alleles  using the upVar and downVar groups




#What proportion of those loci appear in the log2 fold hits?
Varcounts_sig <- as.data.frame(list(sum(length(downVar_hits), length(lowVar_hits), length(highVar_hits), length(upVar_hits)) , length(noVar_hits), length(lowVar_hits), length(highVar_hits), length(upVar_hits), length(downVar_hits)),
                               col.names = c('Var', "noVar", "lowVar", "highVar", "upVar", "downVar"))
Varcounts <- rbind(Varcounts_sig, Varcounts_all[2:7])
for(i in 1:ncol(Varcounts)){
  Varcounts[3,i] <- Varcounts[1,i]/Varcounts[2,i]*100
}
Varcounts
write.xlsx(Varcounts, 'IGR_variation_ratiostest.xlsx')

#Is that difference significant?
#find a test looking at differences in proportions
#chisq test of independence, need a table of format:
#Locus  #Vargroup  #significance
library(reshape2)
Varcounts
contingencytable <- as.data.frame(t(Varcounts[-3,]))

contingencytable$VarGroup <- c('Var', 'NoVar', 'LowVar', 'HighVar', 'UpVar', 'DownVar')
colnames(contingencytable) <- c('sig', 'all', 'VarGroup')
dataset <- melt(contingencytable, id.vars = 'VarGroup')
df <- as.data.frame(matrix(data = NA, nrow = 1723, ncol = 3))
colnames(df) <- c('n', 'VarGroup', 'significance')
df$n <- 1:1723
df$VarGroup[1:244] <- 'Var'
df$VarGroup[245:418] <- 'NoVar'
df$VarGroup[419:1023] <- 'Var'
df$VarGroup[1024:1723] <- 'NoVar'
df$significance[1:418] <- 'significant'
df$significance[419:1723] <- 'nonsignificant'

ggplot(df) +
  aes(x = VarGroup, fill = significance) +
  geom_bar(position = 'fill')

test <- chisq.test(table(df$VarGroup, df$significance))
test
table(df$VarGroup, df$significance)
#What is the average expression for each group?


###Boxplot of hits grouped by high, low, no, up, down var etc on x axis and y axis is the expression. 
#plot the log2 fold change for that locus (between the max and min isolate) for each locus in each group as a box plot or violin plot
#need to make this dataset,
#C1:Locus   C2:log2change   C3:Vargroup  
#nrow = 254
plotdata$log2fold.change <- as.numeric(plotdata$log2fold.change)

ggqqplot(plotdata$log2fold.change[plotdata$VarGroup == 'downVar'])
#Check the data are normally distributed within each group
shapiro.test(plotdata$log2fold.change[plotdata$VarGroup == 'downVar']) #vary the group here to check each
ggqqplot(plotdata$log2fold.change[plotdata$VarGroup == 'lowVar']) #and here
#No groups are normally distributed, so we use non parametric tests to compare means
stats_summary <- group_by(plotdata, VarGroup) %>%
  summarise(
    count = n(),
    median = median(log2fold.change, na.rm = TRUE),
    IQR = IQR(log2fold.change, na.rm = TRUE),
    mean = mean(log2fold.change)
  )

write.xlsx(stats_summary, 'IGR.Groups.Statsummarycorrect.xlsx')

#Below is each pairwise comparison between novar and the other groups
wilcox.test((plotdata$log2fold.change[plotdata$VarGroup == 'downVar']), plotdata$log2fold.change[plotdata$VarGroup == 'NoVar']) #0.077
wilcox.test((plotdata$log2fold.change[plotdata$VarGroup == 'upVar']), plotdata$log2fold.change[plotdata$VarGroup == 'NoVar']) #0.01171
wilcox.test((plotdata$log2fold.change[plotdata$VarGroup == 'highVar']), plotdata$log2fold.change[plotdata$VarGroup == 'NoVar']) #0.3217
wilcox.test((plotdata$log2fold.change[plotdata$VarGroup == 'lowVar']), plotdata$log2fold.change[plotdata$VarGroup == 'NoVar']) #0.0000003014
wilcox.test((plotdata$log2fold.change[plotdata$VarGroup == 'NoMatch']), plotdata$log2fold.change[plotdata$VarGroup == 'NoVar']) #0.5935

wilcox.test((plotdata$log2fold.change[plotdata$VarGroup == 'downVar']), plotdata$log2fold.change[plotdata$VarGroup == 'upVar']) #NS
wilcox.test((plotdata$log2fold.change[plotdata$VarGroup == 'downVar']), plotdata$log2fold.change[plotdata$VarGroup == 'highVar']) #0.0226
wilcox.test((plotdata$log2fold.change[plotdata$VarGroup == 'downVar']), plotdata$log2fold.change[plotdata$VarGroup == 'lowVar']) #NS

wilcox.test((plotdata$log2fold.change[plotdata$VarGroup == 'upVar']), plotdata$log2fold.change[plotdata$VarGroup == 'highVar']) #NS
wilcox.test((plotdata$log2fold.change[plotdata$VarGroup == 'upVar']), plotdata$log2fold.change[plotdata$VarGroup == 'lowVar']) #NS

wilcox.test((plotdata$log2fold.change[plotdata$VarGroup == 'lowVar']), plotdata$log2fold.change[plotdata$VarGroup == 'highVar']) #NS

#Optionally remove 'NoMatch

###
table(plotdata$VarGroup)


mycomparisons <- list(c("downVar", "NoVar"), c("highVar", "NoVar"), c("lowVar", "NoVar"), c("upVar", "NoVar"),  c("highVar", "downVar"), c('lowVar', 'downVar'), c('downVar', 'upVar'), c('upVar', 'highVar'), c('upVar', 'lowVar'))
#Make plot
ggplot(na.omit(plotdata), aes(x = VarGroup, y = log2fold.change, color = VarGroup)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun.y=mean, geom='point', shape = 0, size = 1,  color = 'black') +
  geom_jitter(shape = 16, size = 1.2, position=position_jitter(0.2)) +
  geom_signif(comparisons = mycomparisons, map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05), step_increase = 0.1) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size = 15),
        legend.text = element_text(size=10), legend.title = element_text(size = 15)) +
  labs(y = 'log2-fold change')

#Also need a plot of just Var vs no Var
ggplot(plotdata[!is.na(plotdata$Var),], aes(x = Var, y = log2fold.change, color = Var, na.rm = TRUE)) +
  geom_boxplot(na.rm = T) +
  stat_summary(fun.y=mean, geom='point', shape = 5, size = 2,  color = 'black') +
  geom_jitter(shape = 16, position=position_jitter(0.2)) +
  geom_signif(comparisons = list(c("Var", "NoVar")), map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05)) +
  theme(axis.text = element_text(size=15), axis.title = element_text(size = 20),
        legend.text = element_text(size=15), legend.title = element_text(size = 15))

wilcox.test((plotdata$log2fold.change[plotdata$Var == 'Var']), (plotdata$log2fold.change[plotdata$Var == 'NoVar']))
wilcox.test((plotdata$log2fold.change[plotdata$VarGroup == 'downVar']), plotdata$log2fold.change[plotdata$VarGroup == 'lowVar']) #NS

#####Get the same plots with any logfold change, not just logfold change > 1


###So the two plots show that of the significant RNAseq hits, those that have IGR variance tend to have greater log2 fold changes in expression
#We also observe from the proportion stats that of all loci, those with variation in at least one IGR are overrepresented in significant hits


#Upvar plot
#upplot data structure:
#Locus    Log2foldchange    VarGroup    Allelenumber
alleleplotdata <- alleleplotdata[!is.na(alleleplotdata$Allele.count),]
alleleplotdata$log2fold.change <- as.numeric(alleleplotdata$log2fold.change)
alleleplotdata$Allele.count <- as.numeric(alleleplotdata$Allele.count)

ggplot(alleleplotdata, aes(x=Allele.count, y=log2fold.change)) +
  geom_point(aes(color = VarGroup)) +
  geom_smooth(method = lm, se=F, color = 'black') +
  stat_cor(method = 'pearson') + 
  theme(axis.text = element_text(size=15), axis.title = element_text(size = 20),
        legend.text = element_text(size=15), legend.title = element_text(size = 15))


#Downvar plot


#Combined plot